def complexes(wildcards):
    from glob import glob
    from subprocess import Popen, PIPE
    complexes = []
    for pdbid in glob_wildcards(output_dir + "vorocontacts/{pdbid}.tab").pdbid:
        child = Popen('bin/contact-graph {}vorocontacts/{}.tab --pdb {}pdb/P0DTC2/{}.pdb --output-complexes'.format(output_dir, pdbid, output_dir, pdbid), stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True, text=True, env={'PYTHONPATH': '.'})
        child.stdin.close()
        complexes += [x.strip() for x in child.stdout.readlines()]
    return complexes

rule extract_complexes:
    input:
        complexes
