rule get_skempi_dataset:
    output:
        "dG-datasets/skempi.csv"
    log:
        "dG-datasets/skempi.log"
    shell:
        """
        date --rfc-3339=s > {log}
        wget https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv -O {output} 2>> {log}
        """

rule get_prodigy_dataset:
    output:
        "dG-datasets/prodigy.csv"
    log:
        "dG-datasets/prodigy.log"
    shell:
        """
        date --rfc-3339=s > {log}
        wget https://bianca.science.uu.nl/prodigy/static/PRODIGY_dataset.csv -O {output} 2>> {log}
        """
