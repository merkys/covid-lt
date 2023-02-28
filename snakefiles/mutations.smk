wildcard_constraints:
    chain = "[A-Za-z]",
    mutation = "[A-Z]{2}\d+[A-Z]",
    mutation_maybe_wt = "[A-Z]{2}\d+[A-Z](_wt)?",
    pdbid = "[A-Z0-9]{4}",
    type = "[A-Za-z0-9]+"

def input_complexes():
    complexes = []
    for line in open('SkempiS.txt', 'r').readlines():
        fields = line.split("\t")
        if fields[4] != fields[5]: # Filter out lines where Mutation(s)_PDB and Mutation(s)_cleaned are different
            continue
        if fields[7] != 'forward': # Filter out non-forward (reverse) mutations for now
            continue
        complexes.append("{}_{}{}{}".format(fields[0], fields[4][0], fields[3][0], fields[4][1:]))
        complexes.append("{}_{}{}{}_wt".format(fields[0], fields[4][0], fields[3][0], fields[4][1:]))
    return complexes

rule all_complexes:
    input:
        expand("{complex}.pdb", complex=input_complexes())

rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        mutation = "{pdbid}_{mutation}.pdb",
        wt = "{pdbid}_{mutation}_wt.pdb"
    log:
        "{pdbid}_{mutation}.log"
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
        "{pdbid}_{mutation_maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation_maybe_wt}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            grep ^ATOM {input} \
                | bin/vmd-pdb-to-psf /dev/stdin forcefields/top_all22_prot.rtf \
                | bin/namd-minimize forcefields/par_all22_prot.prm \
                | tar -x --to-stdout output.coor > {output}
        else
            echo -n > {output}
        fi
        """

rule optimize_chain:
    input:
        "{pdbid}_{mutation_maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation_maybe_wt}_{chain}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_select --chain {wildcards.chain} {input} \
            | grep ^ATOM \
            | bin/vmd-pdb-to-psf /dev/stdin forcefields/top_all22_prot.rtf \
            | bin/namd-minimize forcefields/par_all22_prot.prm \
            | tar -x --to-stdout output.coor > {output}
        """

def list_chains(pdb):
    from os import popen
    return popen('grep -h "^ATOM" {} | cut -b 21-22 | uniq'.format(pdb)).read().split()

rule all_optimized:
    input:
        [expand("optimized/{cplx}_{chain}.pdb", cplx=cplx, chain=list_chains(cplx + ".pdb")) for cplx in input_complexes()],
        expand("optimized/{cplx}.pdb", cplx=input_complexes())

rule all_energies:
    input:
        [expand("optimized/{cplx}_{chain}.ener", cplx=cplx, chain=list_chains(cplx + ".pdb")) for cplx in input_complexes()],
        expand("optimized/{cplx}.ener", cplx=input_complexes())
    output:
        solv = "solv.tab",
        vdw = "vdw.tab"
    shell:
        """
        echo -n > {output.solv}
        echo -n > {output.vdw}
        for MUT in $(ls -1 optimized/ | grep -P '^[^_]+_[^_]+\.ener$' | xargs -i basename {{}} .ener | sort | uniq)
        do
            echo -en "$MUT\t" >> {output.solv}
            echo -en "$MUT\t" >> {output.vdw}

            grep -h 'Electrostatic energy' optimized/${{MUT}}.ener optimized/${{MUT}}_wt.ener \
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

            A=$(echo $PARTNERS | cut -c 1)
            B=$(echo $PARTNERS | cut -c 2)

            (
                grep --no-filename '^ENER EXTERN>' optimized/${{MUT}}.ener           | awk '{{print  $3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{MUT}}_${{A}}.ener    | awk '{{print -$3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{MUT}}_${{B}}.ener    | awk '{{print -$3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{MUT}}_wt.ener        | awk '{{print -$3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{MUT}}_wt_${{A}}.ener | awk '{{print  $3}}'
                grep --no-filename '^ENER EXTERN>' optimized/${{MUT}}_wt_${{B}}.ener | awk '{{print  $3}}'
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

rule dssp:
    output:
        part = "sa_part.tab",
        com = "sa_com.tab"
    shell:
        """
        echo -n > {output.part}
        echo -n > {output.com}
        for MUT in $(ls -1 optimized/ | xargs -i basename {{}} .pdb | grep -v _wt | cut -d _ -f 1-2 | cut -d . -f 1 | sort | uniq)
        do
            ORIG_POS=$(echo $MUT | cut -d _ -f 2 | grep -Po '[0-9]+')
            CHAIN=$(echo $MUT | cut -d _ -f 2 | cut -c 2)

            test -s optimized/$(echo $MUT | cut -d _ -f 1)_wt.pdb || continue

            POS=$(grep -P "^${{CHAIN}}${{ORIG_POS}}\s" optimized/$(echo $MUT | cut -d _ -f 1)_wt.map | cut -f 2)

            dssp optimized/$(echo $MUT | cut -d _ -f 1)_wt_$CHAIN.pdb \
                | grep -vP '\.$' \
                | grep -P "^\s+$POS\s" \
                | cut -c 36-38 \
                | xargs -i echo -e $MUT"\t"{{}} >> {output.part} || true
            dssp optimized/$(echo $MUT | cut -d _ -f 1)_wt.pdb \
                | grep -vP '\.$' \
                | grep -P "^\s+$POS\s" \
                | cut -c 36-38 \
                | xargs -i echo -e $MUT"\t"{{}} >> {output.com} || true
        done
        """

rule join_with_skempi:
    input:
        skempi = "SkempiS.txt",
        sa_com = "sa_com.tab",
        sa_part = "sa_part.tab",
        solv = "solv.tab",
        vdw = "vdw.tab"
    output:
        sa_com = "sa_com-skempi.tab",
        sa_part = "sa_part-skempi.tab",
        solv = "solv-skempi.tab",
        vdw = "vdw-skempi.tab"
    shell:
        """
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
