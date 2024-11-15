# Get the container!

You will need to download the PRDM16_CtBP_Project container into
this directory for the workflow scripts to work as designed.

```bash
wget https://BrendelGroup.org/SingularityHub/PRDM16_CtBP_Project.sif
```

Sample usage (our workflow scripts include similar alias lines):

```bash
alias rws="singularity exec -e -B ${PWD}/..  ${PWD}/PRDM16_CtBP_Project.sif"
rws which bedtools
rws which deeptools
rws gth -help
```

Of course this assumes that you have [Apptainer/Singularity](https://apptainer.org/) installed on your system.
Check whether there is package built for your system.
Otherwise, follow the instructions to [install Singularity from source code](https://apptainer.org/user-docs/master/quick_start.html#quick-installation-steps).
