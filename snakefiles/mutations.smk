wildcard_constraints:
    chain = "[A-Za-z]",
    chains = "[A-Za-z]+",
    maybe_wt = "(_wt)?",
    mutation = "[A-Z]{2}-?\d+[a-z]?[A-Z]",
    mutation_maybe_wt = "[A-Z]{2}-?\d+[a-z]?[A-Z](_wt)?",
    pdbid = "[A-Z0-9]{4}",
    type = "[A-Za-z0-9]+"

def skempi_filtered():
    skempi = []
    for line in open('SkempiS.txt', 'r').readlines():
        fields = line.split("\t")
        if fields[7] != 'forward': # Filter out non-forward (reverse) mutations for now
            continue
        skempi.append(fields)
    return skempi

def input_complexes():
    complexes = []
    for fields in skempi_filtered():
        partners = ''.join(sorted(chain[0] for chain in fields[1].split('.') + fields[2].split('.')))
        complexes.append("{}_{}{}{}_{}".format(   fields[0], fields[4][0], fields[3][0], fields[4][1:], partners))
        complexes.append("{}_{}{}{}_{}_wt".format(fields[0], fields[4][0], fields[3][0], fields[4][1:], partners))
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
        "optimized/{pdbid}_{mutation}_{chains}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{chains}{maybe_wt}_{chain}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_select --chain {wildcards.chain} {input} > {output}
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

            test -s optimized/${{STRUCT}}.ener || continue
            test -s optimized/${{STRUCT}}_wt.ener || continue

            echo -en "$MUT\t" >> {output.solv}
            echo -en "$MUT\t" >> {output.vdw}

            grep -h 'Electrostatic energy' optimized/${{STRUCT}}.ener optimized/${{STRUCT}}_wt.ener \
                | awk '{{print $5}}' | xargs echo | awk '{{print $1 - $2}}' >> {output.solv}

            CHAIN=$(echo $MUT | cut -d _ -f 2 | cut -c 2)

            (
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}.ener    | awk '{{print  $3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_wt.ener | awk '{{print -$3}}'
                for CH in $(echo $CHAIN | grep -o .)
                do
                    grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_${{CH}}.ener    | awk '{{print -$3}}'
                    grep --no-filename '^ENER EXTERN>' optimized/${{STRUCT}}_wt_${{CH}}.ener | awk '{{print  $3}}'
                done
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
                | PYTHONPATH=. bin/pdb_charmm_energy /dev/stdin --topology forcefields/top_all36_prot.rtf --parameters forcefields/par_all36m_prot.prm --pbeq \
                | grep -e ^ENER -e 'Electrostatic energy' > {output}
            then
                echo -n > {output}
            fi
        else
            echo -n > {output}
        fi
        """

# FoldX does not understand residues with alternative locations
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
                grep --silent '^DIFF>' $FILE || continue
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

rule N_wt_cont:
    output:
        "N_wt_cont.tab"
    shell:
        """
        for MUT in $(ls -1 optimized/ | grep -P '^[^_]+_[^_]+_[^_]+.pdb$' | xargs -i basename {{}} .pdb | sort | uniq)
        do
            PDB_ID=$(echo $MUT | cut -d _ -f 1)
            CHAIN=$( echo $MUT | cut -d _ -f 2 | cut -c 2)
            CHAINS=$(echo $MUT | cut -d _ -f 3)

            echo -en $MUT"\t"
            bin/pdb_select --chain $CHAINS ${{PDB_ID}}.pdb \
                | PYTHONPATH=. ../bin/distance-contacts /dev/stdin \
                | ../bin/select-contacts --exclude-self --chain $CHAIN \
                | cut -f 1,2,6,7 \
                | awk "{{ if( \$1 == \"$CHAIN\" ) {{print \$3 \$4}} else {{print \$1 \$2}} }}" \
                | sort \
                | uniq \
                | wc -l
        done > {output}
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
        N_wt_cont = "N_wt_cont.tab",
        skempi = "SkempiS.txt"
    output:
        "train-dataset-our.tab"
    shell:
        """
        join {input.vdw} {input.solv} | join - <(sed 's/_[^_]\+\t/\t/' {input.fold}) | join - {input.sa_part} | join - {input.sa_com} | join - <(sed 's/_[^_]\+\t/\t/' {input.N_wt_cont} | sort) | sed 's/ /\t/g' > {output}

        grep forward {input.skempi} \
            | awk '{{print $1 "_" substr($5,1,1) substr($4,1,1) substr($5,2) "\t" $7}}' \
            | sort \
            | join {output} - \
            | cat <(echo mutation vdw solv fold sa_part sa_com cont ddG) - \
            | sed 's/ /\t/g' \
            | sponge {output}
        """

rule train_data_skempi:
    input:
        skempi = "SkempiS.txt"
    output:
        "train-dataset-skempi.tab"
    shell:
        """
        echo mutation vdw solv fold sa_part sa_com cs cont ddG | sed 's/ /\t/g' > {output}

        grep forward {input.skempi} \
            | awk '{{print $1 "_" substr($5,1,1) substr($4,1,1) substr($5,2) "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $7}}' \
            | sed 's/\r//g' \
            | sort >> {output}
        """

rule pdb2pqr:
    input:
        "{name}.pdb"
    output:
        "{name}.pqr"
    singularity:
        "apbs.sif"
    shell:
        "pdb2pqr {input} {output}"

rule apbs:
    input:
        mut = "{name}.pqr",
        wt = "{name}_wt.pqr"
    output:
        "{name}.apbs.out"
    singularity:
        "apbs.sif"
    shell:
        "bin/apbs-pbe {input.mut} {input.wt} > {output} 2>&1"
