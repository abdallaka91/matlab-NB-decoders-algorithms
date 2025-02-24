echo "Bash version ${BASH_VERSION}..."

nb_frame=5
K=21
N=64
GF=64

LANG="en_US.utf8" #define us language to have floating point with point with a point 

echo "start multiple decoding"

for SNR in $(seq -5 1 2)
do
  echo "simulation with SNR= $SNR "
  xterm -xrm '*hold:true' -e ./essai $nb_frame $SNR $K $N $GF&
done 

