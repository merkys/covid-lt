wildcard_constraints:
    chain = "[A-Za-z]",
    chains = "[A-Za-z]+",
    maybe_wt = "(_wt)?",
    mutation = "[A-Z]{2}-?\d+[a-z]?[A-Z]",
    mutation_maybe_wt = "[A-Z]{2}-?\d+[a-z]?[A-Z](_wt)?",
    partner1 = "[A-Za-z]+",
    partner2 = "[A-Za-z]+",
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
        partner1 = ''.join(sorted(chain[0] for chain in fields[1].split('.')))
        partner2 = ''.join(sorted(chain[0] for chain in fields[2].split('.')))
        complexes.append("{}_{}{}{}_{}_{}".format(   fields[0], fields[4][0], fields[3][0], fields[4][1:], partner1, partner2))
        complexes.append("{}_{}{}{}_{}_{}_wt".format(fields[0], fields[4][0], fields[3][0], fields[4][1:], partner1, partner2))
    return complexes

rule all_complexes:
    input:
        expand("{complex}.pdb", complex=input_complexes())

include: "snakefiles/mutations/mutated_complex/FoldX.smk"

rule all_mutated_complex_EvoEF2:
    input:
        expand("EvoEF2/{complex}.pdb", complex=filter(lambda x: not x.endswith("_wt"), input_complexes()))

# EvoEF2 has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex_EvoEF2:
    input:
        "{pdbid}.pdb"
    output:
        "EvoEF2/{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "EvoEF2/{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        bin/pdb_select --first-model {input} \
            | bin/EvoEF2-mutate --mutation {wildcards.mutation} > {output} 2> {log} || echo -n > {output}
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
        "{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    log:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.map"
    shell:
        """
        mkdir --parents $(dirname {output})
        if [ -s {input} ]
        then
            grep ^ATOM {input} \
                | bin/pdb_select --chain {wildcards.partner1}{wildcards.partner2} \
                | PYTHONPATH=. bin/pdb_renumber --output-map {log} \
                | bin/vmd-pdb-to-psf /dev/stdin --topology forcefields/top_all22_prot.rtf --no-split-chains-into-segments 2>/dev/null \
                | bin/namd-minimize forcefields/par_all22_prot.prm 2>/dev/null \
                | tar -x --to-stdout output.coor > {output}
        else
            echo -n > {output}
        fi
        """

def list_chains(name):
    name_parts = name.split('_')
    if name_parts[-1] == 'wt':
        name_parts.pop()
    return list(name_parts[-1])

rule all_optimized:
    input:
        expand("optimized/{cplx}.pdb", cplx=input_complexes())

rule all_charmm_energy:
    input:
        expand("optimized/{cplx}.charmm.ener", cplx=input_complexes())

rule charmm_energy:
    input:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.charmm.ener"
    shell:
        """
        if [ -s {input} ]
        then
            if ! (PYTHONPATH=. bin/pdb_renumber {input} \
                    | PYTHONPATH=. bin/pdb_charmm_energy /dev/stdin --topology forcefields/top_all36_prot.rtf --parameters forcefields/par_all36m_prot.prm --gbsa \
                    | grep -e ^ENER -e '^ Parameter: TOTAL' && \
                  bin/pdb_select --chain {wildcards.partner1} {input} \
                    | PYTHONPATH=. bin/pdb_renumber \
                    | PYTHONPATH=. bin/pdb_charmm_energy /dev/stdin --topology forcefields/top_all36_prot.rtf --parameters forcefields/par_all36m_prot.prm --gbsa \
                    | grep -e ^ENER -e '^ Parameter: TOTAL' && \
                  bin/pdb_select --chain {wildcards.partner2} {input} \
                    | PYTHONPATH=. bin/pdb_renumber \
                    | PYTHONPATH=. bin/pdb_charmm_energy /dev/stdin --topology forcefields/top_all36_prot.rtf --parameters forcefields/par_all36m_prot.prm --gbsa \
                    | grep -e ^ENER -e '^ Parameter: TOTAL') > {output}
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

rule binding_energy_EvoEF:
    input:
        expand("EvoEF2/{complex}.pdb", complex=filter(lambda x: not x.endswith("_wt"), input_complexes()))
    output:
        "binding_energy_EvoEF.tab"
    shell:
        """
        ls -1 EvoEF2/*.pdb \
            | xargs -n 1 basename \
            | while read BASE
              do
                MUTATED=EvoEF2/$BASE
                WT=$(echo $BASE | cut -d _ -f 1).pdb

                test -s $MUTATED || continue
                test -s $WT || continue

                grep --silent ^MODEL $WT && continue || true

                MUT=$(echo $BASE | cut -d _ -f 1-2)

                (
                    EvoEF --command ComputeBinding --split $(echo $BASE | cut -d _ -f 3-4 | sed 's/_/,/g') --pdb $MUTATED
                    EvoEF --command ComputeBinding --split $(echo $BASE | cut -d _ -f 3-4 | sed 's/_/,/g') --pdb $WT
                ) \
                    | grep ^Total \
                    | xargs \
                    | awk '{{print "'$MUT'\t" $3 - $6}}'
              done > {output}
        """

rule binding_energy_EvoEF2:
    input:
        expand("EvoEF2/{complex}.pdb", complex=filter(lambda x: not x.endswith("_wt"), input_complexes()))
    output:
        "binding_energy_EvoEF2.tab"
    shell:
        """
        ls -1 EvoEF2/*.pdb \
            | xargs -n 1 basename \
            | while read BASE
              do
                MUTATED=EvoEF2/$BASE
                WT=$(echo $BASE | cut -d _ -f 1).pdb

                test -s $MUTATED || continue
                test -s $WT || continue

                MUT=$(echo $BASE | cut -d _ -f 1-2)

                (
                    EvoEF2 --command ComputeBinding --split $(echo $BASE | cut -d _ -f 3-4 | sed 's/_/,/g') --pdb $MUTATED
                    EvoEF2 --command ComputeBinding --split $(echo $BASE | cut -d _ -f 3-4 | sed 's/_/,/g') --pdb $WT
                ) \
                    | grep ^Total \
                    | xargs \
                    | awk '{{print "'$MUT'\t" $3 - $6}}'
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

rule apbs_partners:
    input:
        mut = "optimized/{name}_{partner1}_{partner2}.pdb",
        wt = "optimized/{name}_{partner1}_{partner2}_wt.pdb"
    output:
        "optimized/{name}_{partner1}_{partner2}.apbs.solv"
    singularity:
        "apbs.sif"
    shell:
        """
        PYTHONPATH=. bin/apbs-diffeval {wildcards.name} {wildcards.partner1} {wildcards.partner2} > {output}
        """

rule charmm_partners:
    input:
        mut = "optimized/{name}_{partner1}_{partner2}.pdb",
        wt = "optimized/{name_{partner1}_{partner2}_wt.pdb"
    output:
        "optimized/{name}_{partner1}_{partner2}.charmm.solv"
    shell:
        """
        PYTHONPATH=. bin/charmm-diffeval {wildcards.name} {wildcards.partner1} {wildcards.partner2} > {output}
        """

rule all_openmm_energy:
    input:
        expand("optimized/{cplx}.openmm.ener", cplx=input_complexes())

rule existing_openmm_energy:
    output:
        "openmm.tab"
    shell:
        """
        ls -1 optimized/*_wt.openmm.ener \
            | while read FILE
              do
                BASE=$(basename $FILE _wt.openmm.ener)

                echo -en $(echo $BASE | cut -d _ -f 1-2)
                cut -f 1 optimized/${{BASE}}.openmm.ener optimized/${{BASE}}_wt.openmm.ener \
                    | sort \
                    | uniq \
                    | while read FORCE
                      do
                        grep --no-filename ^$FORCE optimized/${{BASE}}.openmm.ener optimized/${{BASE}}_wt.openmm.ener \
                            | cut -f 2 \
                            | xargs \
                            | awk '{{print $1 - $2 - $3 - $4 + $5 + $6}}' \
                            | xargs -i echo -en "\t"{{}}
                      done
                echo
              done > {output}
        """

rule openmm_energy:
    input:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.openmm.ener"
    singularity:
        "container.sif"
    shell:
        """
        (
            sed 's/HSE/HIS/g' {input} \
                | bin/pdb_openmm_minimize --forcefield charmm36.xml --forcefield implicit/gbn2.xml --print-forces --max-iterations 0 --force-unit kcal/mol --split-nonbonded-force
            sed 's/HSE/HIS/g' {input} \
                | bin/pdb_select --chain {wildcards.partner1} \
                | bin/pdb_openmm_minimize --forcefield charmm36.xml --forcefield implicit/gbn2.xml --print-forces --max-iterations 0 --force-unit kcal/mol --split-nonbonded-force
            sed 's/HSE/HIS/g' {input} \
                | bin/pdb_select --chain {wildcards.partner2} \
                | bin/pdb_openmm_minimize --forcefield charmm36.xml --forcefield implicit/gbn2.xml --print-forces --max-iterations 0 --force-unit kcal/mol --split-nonbonded-force
        ) > {output}
        """

rule all_delphi_energy:
    input:
        expand("optimized/{cplx}.delphi.ener", cplx=input_complexes())
    output:
        "delphi.tab"
    shell:
        """
        ls -1 optimized/*_wt.delphi.ener \
            | while read FILE
              do
                BASE=$(basename $FILE _wt.delphi.ener)

                grep --silent '^ Energy> Corrected reaction field energy' optimized/${{BASE}}.delphi.ener || continue
                grep --silent '^ Energy> Corrected reaction field energy' optimized/${{BASE}}_wt.delphi.ener || continue

                DIFF=$(grep --no-filename '^ Energy> Corrected reaction field energy' optimized/${{BASE}}.delphi.ener optimized/${{BASE}}_wt.delphi.ener \
                    | awk '{{print $7}}' \
                    | xargs \
                    | awk '{{print $1 - $2 - $3 - $4 + $5 + $6}}')
                echo -e $(echo $BASE | cut -d _ -f 1-2)"\t"${{DIFF}}
              done | tee {output}
        """

rule delphi_energy:
    input:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.delphi.ener"
    shell:
        """
        (
            bin/pdb_delphi_energy {input}
            bin/pdb_select --chain {wildcards.partner1} {input} | bin/pdb_delphi_energy
            bin/pdb_select --chain {wildcards.partner2} {input} | bin/pdb_delphi_energy
        ) > {output}
        """

rule all_sander_energy:
    input:
        expand("optimized/{cplx}.sander.ener", cplx=input_complexes())

rule sander_energy:
    input:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}{maybe_wt}.sander.ener"
    shell:
        """
        (
            grep ^ATOM {input} | bin/tleap --output-tar --source leaprc.protein.ff19SB 2>/dev/null | bin/sander --input-tar
            grep ^ATOM {input} | bin/pdb_select --chain {wildcards.partner1} | bin/tleap --output-tar --source leaprc.protein.ff19SB 2>/dev/null | bin/sander --input-tar
            grep ^ATOM {input} | bin/pdb_select --chain {wildcards.partner2} | bin/tleap --output-tar --source leaprc.protein.ff19SB 2>/dev/null | bin/sander --input-tar
        ) > {output}
        """

rule all_UEP:
    input:
        expand("optimized/{complex}.uep.csv", complex=filter(lambda x: x.endswith("_wt"), input_complexes()))

rule existing_UEP:
    output:
        "uep.tab"
    shell:
        """
        ls -1 optimized/*.uep.csv \
            | while read FILE
              do
                PDB=$(basename $FILE | cut -d _ -f 1)
                MUT=$(basename $FILE | cut -d _ -f 2)
                DDG=$(bin/process-uep $FILE --map optimized/$(basename $FILE .uep.csv).map --mutation $MUT)
                echo $DDG | grep --silent . && echo -e ${{PDB}}_${{MUT}}"\t"$DDG || true
              done | tee {output}
        """

rule UEP:
    input:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    output:
        "optimized/{pdbid}_{mutation}_{partner1}_{partner2}_wt.uep.csv"
    singularity:
        "UEP.sif"
    shell:
        """
        TMPFILE=$(mktemp --suffix .pdb)
        sed 's/HSE/HIS/g' {input} > $TMPFILE
        PYTHONPATH=UEP python3 UEP/UEP.py --pdb $TMPFILE --interface {wildcards.partner1},{wildcards.partner2}
        rm -f $TMPFILE
        mv /tmp/$(basename $TMPFILE .pdb)_UEP_*.csv {output}
        """
