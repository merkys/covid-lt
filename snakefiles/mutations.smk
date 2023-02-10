wildcard_constraints:
    chain = "[A-Za-z]",
    mutation = "[A-Z]{2}\d+[A-Z]",
    pdbid = "[A-Z0-9]{4}",
    type = "[A-Za-z0-9]+"

def input_complexes():
    complexes = []
    seen_wt = set()
    for line in open('SkempiS.txt', 'r').readlines():
        fields = line.split("\t")
        if fields[4] != fields[5]: # Filter out lines where Mutation(s)_PDB and Mutation(s)_cleaned are different
            continue
        if fields[7] != 'forward': # Filter out non-forward (reverse) mutations for now
            continue
        if fields[0] not in seen_wt:
            complexes.append("{}_wt".format(fields[0]))
            seen_wt.add(fields[0])
        complexes.append("{}_{}{}{}".format(fields[0], fields[4][0], fields[3][0], fields[4][1:]))
    return complexes

rule all_complexes:
    input:
        expand("{complex}.pdb", complex=input_complexes())

rule mutated_complex:
    input:
        "{pdbid}_wt.pdb"
    output:
        "{pdbid}_{mutation}.pdb"
    shell:
        """
        if [ -s {input} ]
        then
            bin/EvoEF2-mutate {input} {wildcards.mutation} > {output} || echo -n > {output}
        else
            echo -n > {output}
        fi
        """

rule wt:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_wt.pdb"
    shell:
        """
        PYTHONPATH=. bin/pdb_clean_incomplete {input} > {output}
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
        "{pdbid}_{type}.pdb"
    output:
        "optimized/{pdbid}_{type}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            pdb_renumber --from 1 {input} \
                | bin/vmd-pdb-to-psf /dev/stdin forcefields/top_all22_prot.rtf \
                | bin/namd-minimize forcefields/par_all22_prot.prm \
                | tar -x --to-stdout output.coor > {output}
        else
            echo -n > {output}
        fi
        """

rule optimize_chain:
    input:
        "{pdbid}_{type}.pdb"
    output:
        "optimized/{pdbid}_{type}_{chain}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_select --chain {wildcards.chain} {input} \
            | pdb_renumber --from 1 \
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
        for MUT in $(ls -1 optimized/ | xargs -i basename {{}} .ener | cut -d _ -f 1-2 | cut -d . -f 1 | sort | uniq)
        do
            test -s optimized/$MUT.ener || continue
            test -s optimized/$(echo $MUT | cut -d _ -f 1)_wt.ener || continue
            test $(echo $MUT | cut -d _ -f 2) == wt && continue

            for OUT in {output.solv} {output.vdw}
            do
                echo -n $(echo $MUT | sed 's/_/\t/g') >> $OUT
            done

            (
                grep 'Electrostatic energy' optimized/$(echo $MUT | cut -d _ -f 1)_wt.ener | awk '{{print $NF}}'
                grep --no-filename 'Electrostatic energy' optimized/$(echo $MUT | cut -d _ -f 1)_wt_*.ener | awk '{{print -$NF}}'
            ) | bin/sum | xargs -i echo -n "\t"{{}} >> {output.solv}

            (
                grep '^ENER EXTERN>' optimized/$(echo $MUT | cut -d _ -f 1)_wt.ener | awk '{{print $3}}'
                grep --no-filename '^ENER EXTERN>' optimized/$(echo $MUT | cut -d _ -f 1)_wt_*.ener | awk '{{print -$3}}'
            ) | bin/sum | xargs -i echo -n "\t"{{}} >> {output.vdw}

            echo -n "\t"$(grep 'Electrostatic energy' optimized/$MUT.ener | awk '{{print $NF}}') >> {output.solv}
            echo -n "\t"$(grep '^ENER EXTERN>' optimized/$MUT.ener | awk '{{print $3}}') >> {output.vdw}

            for CHAIN in optimized/${{MUT}}_*.ener
            do
                echo -n "\t"$(grep 'Electrostatic energy' $CHAIN | awk '{{print $NF}}') >> {output.solv}
                echo -n "\t"$(grep '^ENER EXTERN>' $CHAIN | awk '{{print $3}}') >> {output.vdw}
            done
            echo >> {output.solv}
            echo >> {output.vdw}
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
            if ! pdb_renumber --from 1 {input} \
                | bin/pdb_charmm_energy /dev/stdin --topology forcefields/top_all22_prot.rtf --parameters forcefields/par_all22_prot.prm --pbeq \
                | grep -e ^ENER -e 'Electrostatic energy' > {output}
            then
                echo -n > {output}
            fi
        else
            echo -n > {output}
        fi
        """

rule join_with_skempi:
    input:
        skempi = "SkempiS.txt",
        solv = "solv.tab",
        vdw = "vdw.tab"
    output:
        solv = "solv-skempi.tab",
        vdw = "vdw-skempi.tab"
    shell:
        """
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $14}} }}' | sort -k1.1) \
             <(sed 's/ /_/' {input.solv} | sort -k1.1) | sed 's/ /\t/g' > {output.solv}
        join <(tail -n +2 {input.skempi} | awk '{{ if( $8 == "forward" ) {{print $1 "_" substr($5,0,1) substr($4,0,1) substr($5,2) "\t" $13}} }}' | sort -k1.1) \
             <(sed 's/ /_/' {input.vdw}  | sort -k1.1) | sed 's/ /\t/g' > {output.vdw}
        """
