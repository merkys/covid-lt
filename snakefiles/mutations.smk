wildcard_constraints:
    chain = "[A-Za-a]",
    mutation = "[A-Za-z0-9]+"

def chains():
    from os import popen
    return popen('grep -hE "^ATOM|^HETATM" wt.pdb | cut -b 21-22 | uniq').read().split()

def mutations():
    return [x.strip() for x in open('mutations.lst', 'r').readlines()]

checkpoint generate_mutations:
    input:
        wt = 'wt.pdb',
        mutations = 'mutations.lst'
    output:
        expand("{mut}.pdb", mut=mutations())
    shell:
        """
        EvoEF2 --command BuildMutant --pdb {input.wt} --mutant_file {input.mutations}
        paste <(ls -1 wt_Model_00*) <(echo {output} | xargs -n1 echo) | while read LINE; do mv $LINE; done
        """

rule energies:
    input:
        expand('complex/{mut}.ener', mut=mutations()),
        expand('{chain}/{mut}.ener', chain=chains(), mut=mutations()),
        'complex/wt.ener',
        expand('{chain}/wt.ener', chain=chains())
    output:
        "vdw.tab"
    shell:
        """
        echo -n > {output}
        for MUT in $(ls -1 {input} | xargs -i basename {{}} .ener | sort | uniq)
        do
            echo -n $MUT >> {output}
            for CHAIN in complex $(find . -name '[A-Z]' | sort)
            do
                echo -n "\t"$(grep '^ENER EXTERN>' $CHAIN/$MUT.ener | awk '{{print $3}}') >> {output}
            done
            echo >> {output}
        done
        """

rule chain:
    input:
        '{mutation}.pdb'
    output:
        '{chain}/{mutation}.pdb'
    shell:
        """
        mkdir --parents $(dirname {output})
        ../bin/pdb_select --chain {wildcards.chain} {input} | pdb_renumber --from 1 | ../bin/vmd-pdb-to-psf /dev/stdin ../1A22.namd/top_all22_prot.rtf | ../bin/namd-minimize ../1A22.namd/par_all22_prot.prm | tar -x --to-stdout output.coor > {output}
        """

rule complex:
    input:
        '{mutation}.pdb'
    output:
        'complex/{mutation}.pdb'
    shell:
        """
        mkdir --parents $(dirname {output})
        ../bin/vmd-pdb-to-psf {input} ../1A22.namd/top_all22_prot.rtf | ../bin/namd-minimize ../1A22.namd/par_all22_prot.prm | tar -x --to-stdout output.coor | pdb_renumber --from 1 > {output}
        """

rule energy:
    input:
        '{name}.pdb'
    output:
        '{name}.ener'
    shell:
        "../bin/pdb_charmm_energy {input} ../1A22.namd/top_all22_prot.rtf ../1A22.namd/par_all22_prot.prm | grep ^ENER > {output}"
