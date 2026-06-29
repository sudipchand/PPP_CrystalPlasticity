#!/bin/bash
#PBS -N cp_umat
#PBS -q teachingq
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -l walltime=04:00:00
#PBS -o cp_umat.out
#PBS -e cp_umat.err

echo "Job started from $(pwd)"
cd $PBS_O_WORKDIR
echo "Working directory: $(pwd)"

# ---------------------------------------------------------------
# Load modules (same as working elastic UMAT example)
# ---------------------------------------------------------------
module purge
module load gcc/11.4.0
module load intel/2024.2.0
module load abaqus/2023

echo "Loaded modules:"
module list
echo "Paths:"
which ifort
which abaqus

# ---------------------------------------------------------------
# Job variables
# ---------------------------------------------------------------
JOBNAME=cp_uniaxial
UMATFILE=UMAT_crystal_plasticity_aniso_elast_interact_hard_wo_module.f
NCPUS=4   # must match mpiprocs above

# ---------------------------------------------------------------
# Clean up old output files from previous runs
# ---------------------------------------------------------------
rm -f ${JOBNAME}.com ${JOBNAME}.dat ${JOBNAME}.msg ${JOBNAME}.odb
rm -f ${JOBNAME}.prt ${JOBNAME}.sim ${JOBNAME}.sta ${JOBNAME}.log
rm -f ${JOBNAME}.lck ${JOBNAME}.feh

# ---------------------------------------------------------------
# Run Abaqus
# ---------------------------------------------------------------
abaqus job=${JOBNAME} user=${UMATFILE} cpus=${NCPUS} interactive

echo "Abaqus finished."

# Quick post-run checks
if [ -f ${JOBNAME}.sta ]; then
    grep -i "completed\|error\|aborted" ${JOBNAME}.sta
fi
if [ -f ${JOBNAME}.dat ]; then
    grep -i "error" ${JOBNAME}.dat || echo "No errors in ${JOBNAME}.dat"
fi
if [ -f ${JOBNAME}.feh ]; then
    echo "--- UMAT log (${JOBNAME}.feh) ---"
    cat ${JOBNAME}.feh
fi

echo "Job ended."
