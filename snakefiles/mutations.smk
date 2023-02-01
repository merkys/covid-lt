wildcard_constraints:
    mutation = "[A-Za-z0-9]+"

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

rule complex_energies:
    input:
        expand('complex/{mut}.ener', mut=mutations()),
        'complex/wt.ener'

rule complex:
    input:
        '{mutation}.pdb'
    output:
        'complex/{mutation}.pdb'
    shell:
        """
        mkdir --parents $(dirname {output})
        ../bin/vmd-pdb-to-psf {input} ../1A22.namd/top_all22_prot.rtf | ../bin/namd-minimize ../1A22.namd/par_all22_prot.prm | tar -x --to-stdout output.coor > {output}
        """

rule energy:
    input:
        '{name}.pdb'
    output:
        '{name}.ener'
    shell:
        "../bin/pdb_charmm_energy {input} ../1A22.namd/top_all22_prot.rtf ../1A22.namd/par_all22_prot.prm | grep ^ENER > {output}"
