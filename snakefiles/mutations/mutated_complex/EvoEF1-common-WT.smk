# EvoEF has trouble dealing with multi-model PDB files: https://github.com/xiaoqiah/EvoEF2/issues/2
rule mutated_complex:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}.pdb"
    log:
        "{pdbid}_{mutation}_{partner1}_{partner2}.log"
    shell:
        """
        if ! test -s {input}
        then
            echo -n > {output}
            exit
        fi

        if ! bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input} \
            | bin/EvoEF2-mutate --mutation {wildcards.mutation} --EvoEF2-command EvoEF > {output} 2> {log}
        then
            echo -n > {output}
        fi
        """

rule wild_type:
    input:
        "{pdbid}.pdb"
    output:
        "{pdbid}_{mutation}_{partner1}_{partner2}_wt.pdb"
    shell:
        """
        if ls -1 {wildcards.pdbid}_*_{wildcards.partner1}_{wildcards.partner2}_wt.pdb | grep --silent .
        then
            # Already generated, sufficient to copy
            ls -1 {wildcards.pdbid}_*_{wildcards.partner1}_{wildcards.partner2}_wt.pdb \
                | head -n 1 \
                | xargs -i cp {{}} {output}
        else
            voronota-contacts -i {input} \
                | bin/vorocontacts2tab \
                | bin/select-contacts --between {wildcards.partner1} --between {wildcards.partner2} --exclude-self --exclude-hetero \
                | cut -f 1-3,6-8 \
                | while read CHAIN1 POS1 RESIDUE1 CHAIN2 POS2 RESIDUE2
                  do
                    echo $(echo $RESIDUE1 | bin/aa-3to1)$CHAIN1$POS1$(echo $RESIDUE1 | bin/aa-3to1)
                    echo $(echo $RESIDUE2 | bin/aa-3to1)$CHAIN2$POS2$(echo $RESIDUE2 | bin/aa-3to1)
                  done \
                | sort \
                | uniq \
                | xargs echo \
                | sed 's/ /,/g' \
                | xargs -i bin/EvoEF2-mutate <(bin/pdb_select --first-model --chain {wildcards.partner1}{wildcards.partner2} {input}) \
                    --mutation {{}} --EvoEF2-command EvoEF --wt {output} > /dev/null
        fi
        """
