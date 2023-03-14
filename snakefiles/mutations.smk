wildcard_constraints:
    chain = "[A-Za-z]",
    chains = "[A-Za-z]+",
    maybe_wt = "(_wt)?",
    mutation = "[A-Z]{2}\d+[A-Z]",
    mutation_maybe_wt = "[A-Z]{2}\d+[A-Z](_wt)?",
    pdbid = "[A-Z0-9]{4}",
    type = "[A-Za-z0-9]+"

def skempi_filtered():
    skempi = []
    for line in open('SkempiS.txt', 'r').readlines():
        fields = line.split("\t")
        if fields[1].count('.') > 0 or fields[2].count('.') > 0: # Filter out lines having more than two partners
            continue
        if fields[4] != fields[5]: # Filter out lines where Mutation(s)_PDB and Mutation(s)_cleaned are different
            continue
        if fields[7] != 'forward': # Filter out non-forward (reverse) mutations for now
            continue
        skempi.append(fields)
    return skempi

def input_complexes():
    complexes = []
    for fields in skempi_filtered():
        complexes.append("{}_{}{}{}_{}{}".format(   fields[0], fields[4][0], fields[3][0], fields[4][1:], fields[1][0], fields[2][0]))
        complexes.append("{}_{}{}{}_{}{}_wt".format(fields[0], fields[4][0], fields[3][0], fields[4][1:], fields[1][0], fields[2][0]))
    return complexes

def skempi_get_details(mutation):
    if mutation[1].isalpha(): # chain given
        mutation = mutation[0] + mutation[2:]
    for fields in skempi_filtered():
        if fields[4] == mutation:
            return fields
    return None

rule all_complexes:
    input:
        expand("{complex}.pdb", complex=input_complexes())

rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        mutation = "{pdbid}_{mutation}_{chains}.pdb",
        wt = "{pdbid}_{mutation}_{chains}_wt.pdb"
    log:
        "{pdbid}_{mutation}_{chains}.log"
    singularity:
        "container.sif"
    shell:
        """
        echo -n > {output.mutation}
        echo -n > {output.wt}
        if [ -s {input} ]
        then
            PATH=~/bin:$PATH bin/FoldX-mutate {input} {wildcards.mutation} rotabase.txt {output.wt} > {output.mutation} 2> {log} || echo Processing {input} failed >&2
        fi
        """

rule all_original_pdbs:
    input:
        expand("{id}.pdb", id=list(set([x[0:4] for x in input_complexes()])))

rule original_pdb:
    output:
        "{pdbid}.pdb"
    shell:
        """
        wget https://files.rcsb.org/download/{output} -O {output} || echo PDB file for {output} cannot be downloaded >&2
        chmod -w {output} 2>/dev/null || true # Intentional
        """

rule optimize_complex:
    input:
        "{pdbid}_{mutation}_{chains}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{chains}{maybe_wt}.pdb"
    log:
        "optimized/{pdbid}_{mutation}_{chains}{maybe_wt}.map"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            grep ^ATOM {input} \
                | bin/pdb_select --chain {wildcards.chains} \
                | PYTHONPATH=. bin/pdb_renumber --output-map {log} \
                | bin/vmd-pdb-to-psf /dev/stdin --topology forcefields/top_all22_prot.rtf --no-split-chains-into-segments \
                | bin/namd-minimize forcefields/par_all22_prot.prm \
                | tar -x --to-stdout output.coor > {output}
        else
            echo -n > {output}
        fi
        """

rule optimize_chain:
    input:
        "{pdbid}_{mutation}_{chains}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{chains}{maybe_wt}_{chain}.pdb"
    log:
        "optimized/{pdbid}_{mutation}_{chains}{maybe_wt}_{chain}.map"
    shell:
        """
        mkdir --parents $(dirname {output})
        grep ^ATOM {input} \
            | bin/pdb_select --chain {wildcards.chain} \
            | PYTHONPATH=. bin/pdb_renumber --output-map {log} \
            | bin/vmd-pdb-to-psf /dev/stdin --topology forcefields/top_all22_prot.rtf --no-split-chains-into-segments \
            | bin/namd-minimize forcefields/par_all22_prot.prm \
            | tar -x --to-stdout output.coor > {output}
        """

def list_chains(name):
    name_parts = name.split('_')
    if name_parts[-1] == 'wt':
        name_parts.pop()
    return list(name_parts[-1])

rule all_optimized:
    input:
        [expand("optimized/{cplx}_{chain}.pdb", cplx=cplx, chain=list_chains(cplx)) for cplx in input_complexes()],
        expand("optimized/{cplx}.pdb", cplx=input_complexes())

rule all_energies:
    input:
        [expand("optimized/{cplx}_{chain}.ener", cplx=cplx, chain=list_chains(cplx)) for cplx in input_complexes()],
        expand("optimized/{cplx}.ener", cplx=input_complexes())
    output:
        solv = "solv.tab",
        vdw = "vdw.tab"
    shell:
        """
        echo -n > {output.solv}
        echo -n > {output.vdw}
        for STRUCT in $(ls -1 optimized/ | grep -P '^[^_]+_[^_]+_[^_]+\.ener$' | xargs -i basename {{}} .ener | sort | uniq)
        do
            MUT=$(echo $STRUCT | cut -d _ -f 1-2)

            echo -en "$MUT\t" >> {output.solv}
            echo -en "$MUT\t" >> {output.vdw}

            grep -h 'Electrostatic energy' optimized/${{STRUCT}}.ener optimized/${{STRUCT}}_wt.ener \
                | awk '{{print $5}}' | xargs echo | awk '{{print $1 - $2}}' >> {output.solv}

            PDB=$(echo $MUT | cut -d _ -f 1)
            CHAIN=$(echo $MUT | cut -d _ -f 2 | cut -c 2)
            CHAINLESS_MUT=$(echo $MUT | cut -d _ -f 2 | cut -c 1,3-)

            PARTNERS=$(grep ^$PDB SkempiS.txt \
                        | grep forward \
                        | awk '{{if( $5 == $6 )                 {{print $0}}}}' \
                        | awk '{{if( $4 == "'$CHAIN'_1" )       {{print $0}}}}' \
                        | awk '{{if( $5 == "'$CHAINLESS_MUT'" ) {{print $0}}}}' \
                        | awk '{{print substr($2,1,1) substr($3,1,1)}}' \
                        | head -n1)

            A=$(echo $STRUCT | cut -d _ -f 3 | cut -c 1)
            B=$(echo $STRUCT | cut -d _ -f 3 | cut -c 2)

            (
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}.ener           | awk '{{print  $3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_${{A}}.ener    | awk '{{print -$3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_${{B}}.ener    | awk '{{print -$3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_wt.ener        | awk '{{print -$3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_wt_${{A}}.ener | awk '{{print  $3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_wt_${{B}}.ener | awk '{{print  $3}}'
            ) | bin/sum >> {output.vdw}
        done
        """

rule energy:
    input:
        "{name}.pdb"
    output:
        "{name}.ener"
    shell:
        """
        if [ -s {input} ]
        then
            if ! PYTHONPATH=. bin/pdb_renumber {input} \
                | bin/pdb_charmm_energy /dev/stdin --topology forcefields/top_all22_prot.rtf --parameters forcefields/par_all22_prot.prm --pbeq \
                | grep -e ^ENER -e 'Electrostatic energy' > {output}
            then
                echo -n > {output}
            fi
        else
            echo -n > {output}
        fi
        """

rule fold_energy:
    input:
        expand("{complex}.log", complex=filter(lambda x: not x.endswith("_wt"), input_complexes()))
    output:
        "fold.tab"
    shell:
        """
        ls -1 *.log \
            | while read FILE
              do
                echo -en $(basename $FILE .log)"\t"
                grep '^DIFF>' $FILE | tail -n 1 | cut -f 2
              done > {output}
        """

rule dssp:
    output:
        part = "sa_part.tab",
        com = "sa_com.tab"
    shell:
        """
        echo -n > {output.part}
        echo -n > {output.com}
        for MUT in $(ls -1 optimized/ | grep -P '^[^_]+_[^_]+_[^_]+.pdb$' | xargs -i basename {{}} .pdb | sort | uniq)
        do
            ORIG_POS=$(echo $MUT | cut -d _ -f 2 | grep -Po '[0-9]+')
            CHAIN=$(echo $MUT | cut -d _ -f 2 | cut -c 2)

            test -s optimized/${{MUT}}_wt.pdb || continue

            POS=$(grep -P "^${{CHAIN}}${{ORIG_POS}}\s" optimized/${{MUT}}_wt.map | cut -f 2)
            POS_IN_CHAIN=$(grep -P "^${{CHAIN}}${{ORIG_POS}}\s" optimized/${{MUT}}_wt_$CHAIN.map | cut -f 2)

            dssp optimized/${{MUT}}_wt_$CHAIN.pdb \
                | grep -vP '\.$' \
                | grep -P "^\s+$POS_IN_CHAIN\s" \
                | cut -c 36-38 \
                | xargs -i echo -e $(echo $MUT | cut -d _ -f 1-2)"\t"{{}} >> {output.part} || true
            dssp optimized/${{MUT}}_wt.pdb \
                | grep -vP '\.$' \
                | grep -P "^\s+$POS\s" \
                | cut -c 36-38 \
                | xargs -i echo -e $(echo $MUT | cut -d _ -f 1-2)"\t"{{}} >> {output.com} || true
        done
        """

rule join_with_skempi:
    input:
        skempi = "SkempiS.txt",
        fold = "fold.tab",
        sa_com = "sa_com.tab",
        sa_part = "sa_part.tab",
        solv = "solv.tab",
        vdw = "vdw.tab"
    output:
        fold = "fold-skempi.tab",
        sa_com = "sa_com-skempi.tab",
        sa_part = "sa_part-skempi.tab",
        solv = "solv-skempi.tab",
        vdw = "vdw-skempi.tab"
    shell:
        """
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $15}} }}' | sort -k1.1) \
             <(sort -k1.1 {input.fold}) | sed 's/ /\t/g' > {output.fold}
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $16}} }}' | sort -k1.1) \
             <(sort -k1.1 {input.sa_com}) | sed 's/ /\t/g' > {output.sa_com}
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $17}} }}' | sort -k1.1) \
             <(sort -k1.1 {input.sa_part}) | sed 's/ /\t/g' > {output.sa_part}
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $14}} }}' | sort -k1.1) \
             <(sed 's/ /_/' {input.solv} | sort -k1.1) | sed 's/ /\t/g' > {output.solv}
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $13}} }}' | sort -k1.1) \
             <(sed 's/ /_/' {input.vdw}  | sort -k1.1) | sed 's/ /\t/g' > {output.vdw}
        """

rule side_by_side:
    input:
        "{type}-skempi.tab"
    output:
        "{type}-diff.tab"
    shell:
        """
        cat {input} \
            | while read LINE
                do
                    echo -n $(echo "$LINE" | cut -f 1-2)"\t"
                    Rscript -e "0 + $(echo "$LINE" | cut -f 3) - ($(echo "$LINE" | cut -f 4) - $(echo "$LINE" | cut -f 5- | ./bin/sum))" | cut -d ' ' -f 2
                done \
            | sed 's/ /\t/g' > {output}
        """

rule train_dataset_our:
    input:
        vdw = "vdw.tab",
        solv = "solv.tab",
        fold = "fold.tab",
        sa_part = "sa_part.tab",
        sa_com = "sa_com.tab",
        skempi = "SkempiS.txt"
    output:
        "train-dataset-our.tab"
    shell:
        """
        join {input.vdw} {input.solv} | join - <(sed 's/_AB\t/\t/' {input.fold}) | join - {input.sa_part} | join - {input.sa_com} | sed 's/ /\t/g' > {output}

        grep forward {input.skempi} \
            | awk '{{if( $5 == $6 )   {{print $0}}}}' \
            | awk '{{if( $2 !~ /\./ ) {{print $0}}}}' \
            | awk '{{if( $3 !~ /\./ ) {{print $0}}}}' \
            | awk '{{print $1 "_" substr($5,1,1) substr($4,1,1) substr($5,2) "\t" $7}}' \
            | sort \
            | join {output} - \
            | cat <(echo mutation vdw solv fold sa_part sa_com ddG) - \
            | sed 's/ /\t/g' \
            | sponge {output}
        """
