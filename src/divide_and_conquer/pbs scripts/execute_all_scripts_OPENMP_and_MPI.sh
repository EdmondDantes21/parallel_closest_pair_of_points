for file in ./MPI_and_OPENMP/*
do
  qsub "$file"
done