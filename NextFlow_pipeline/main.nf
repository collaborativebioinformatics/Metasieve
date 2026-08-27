#!/usr/bin/env nextflow

process RunPythonScript {
    input:
    path my_script

    output:
    stdout

    script:
    """
    python3 ${my_script}
    """
}

workflow {
    script_ch = file('script.py')
    RunPythonScript(script_ch) | view
}
