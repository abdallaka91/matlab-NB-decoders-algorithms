echo "Bash version ${BASH_VERSION}..."

nb_frame=9999999999

K=42
N=128
GF=64
nm=20
offset=0.6

LANG="en_US.utf8" #define us language to have floating point with point with a point 

echo "start multiple decoding"

for SNR in $(seq -1.6 0.4 0.8)
do
  echo "simulation with SNR= $SNR "
  xterm -xrm '*hold:true' -e ./essai $nb_frame $SNR $K $N $GF $nm $offset&
done 

