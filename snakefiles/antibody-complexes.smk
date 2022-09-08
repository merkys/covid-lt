def complexes():
    from glob import glob
    from subprocess import Popen, PIPE
    complexes = []
    for pdbid in glob_wildcards(output_dir + "vorocontacts/{pdbid}.tab").pdbid:
        child = Popen('bin/contact-graph {}vorocontacts/{}.tab --pdb {}pdb/P0DTC2/{}.pdb --output-complexes'.format(output_dir, pdbid, output_dir, pdbid), stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True, text=True, env={'PYTHONPATH': '.'})
        child.stdin.close()
        complexes += [pdbid + '_' + x.strip() for x in child.stdout.readlines()]
    return complexes

rule extract_complexes:
    input:
        expand(output_dir + "pdb/antibodies/complexes/{complex}.pdb", complex=complexes())
    output:
        output_dir + "pdb/antibodies/complexes/out.log"
    shell:
        "touch {output}"

rule extract_complex:
    input:
        output_dir + "pdb/P0DTC2/{pdbid}.pdb"
    output:
        output_dir + "pdb/antibodies/complexes/{pdbid}_{chains}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        bin/pdb_select --chain $(echo {wildcards.chains} | cut -c 1) --chain $(echo {wildcards.chains} | cut -c 2) --chain $(echo {wildcards.chains} | cut -c 3) {input} \
            | PYTHONPATH=. bin/pdb_rename_chains --map $(echo {wildcards.chains} | cut -c 1):A --map $(echo {wildcards.chains} | cut -c 2):H --map $(echo {wildcards.chains} | cut -c 3):L > {output}
        """

rule renumbered:
    input:
        output_dir + "pdb/antibodies/complexes/{name}.pdb"
    output:
        output_dir + "pdb/antibodies/renumbered/{name}.pdb"
    shell:
        """
        mkdir --parents $(dirname {output})
        convert_pdb_to_antibody_numbering_scheme.py {input} {output} H L c
        """

rule snugdock:
    input:
        output_dir + "pdb/antibodies/renumbered/{name}.pdb"
    output:
        output_dir + "pdb/antibodies/snugdock/{name}.pdb"
    log:
        output_dir + "pdb/antibodies/snugdock/{name}.log"
    shell:
        """
        mkdir --parents $(dirname {output})
        TMP_DIR=$(mktemp --directory)
        (cd $TMP_DIR && snugdock -s {input} -partners LH_A) > {log}
        mv $TMP_DIR/{wildcards.name}_0001.pdb {output}
        rm -rf $TMP_DIR
        """
