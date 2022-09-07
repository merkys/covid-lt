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
        expand(output_dir + "pdb/split/{complex}.pdb", complex=complexes())

rule extract_complex:
    input:
        output_dir + "pdb/P0DTC2/{pdbid}.pdb"
    output:
        output_dir + "pdb/split/{pdbid}_{chains}.pdb"
    shell:
        """
        bin/pdb_select --chain $(echo {wildcards.chains} | cut -c 1) $(echo {wildcards.chains} | cut -c 2) $(echo {wildcards.chains} | cut -c 3) {input} \
            | PYTHONPATH=. bin/pdb_rename_chains --map $(echo {wildcards.chains} | cut -c 1):A --map $(echo {wildcards.chains} | cut -c 2):H --map $(echo {wildcards.chains} | cut -c 3):L > {output}
        """
