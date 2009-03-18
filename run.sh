#!/bin/bash -e

echo;echo;echo "Hello $USER";echo 'Welcome to the automatic mesh generator';echo

#inimesh='datafiles/mesh'
#plot=T
#ang=20
read inimesh plot ang < "specsfile"

echo Making Components...
OMP_NUM_THREADS=2
make

echo;echo 'Generating initial mesh'
cp domain meshgen/
cp holefile meshgen/
cd meshgen
    ./meshgen mesh.poly
    mv mesh.poly ../datafiles/
    echo Created $inimesh.poly;echo
    cp domain.out ../adaptive/domain
cd ..
    
cp domain nssim/domain
sed '1d' holefile > adaptive/holefile

echo Calling triangle to triangulate region...
triangle/./triangle -co2Qq$ang $inimesh.poly
if [ $plot == 'T' ] ; then 
    triangle/./showme $inimesh.1.ele 
fi

echo Created $inimesh.1.node
echo Created $inimesh.1.poly
echo Created $inimesh.1.ele;echo

echo Formatting triangle output...
sed '1d' $inimesh.1.ele > nssim/elements.in
sed '1d' $inimesh.1.node > nssim/nodes.in
echo

echo '1' > summary.out
cp simfile nssim

cd nssim
	echo;echo Running 1st NS Simulation..
	#shark -i -G -T 1s ./nssim  #Use this line for code profiling
	./nssim elements.in nodes.in
	cat summary.out >> ../summary.out
    cat time.out >> ../time.out
cd ..;echo

echo ok > adaptive/continue
cp refinefile adaptive/
cd adaptive
    ./adapt ../nssim/nodes.out ../nssim/elements.out  \
    ../nssim/solution.out ../nssim/corner.out mesh.1.poly \
    read.table.o
    cat summary.out >> ../summary.out
cd ..

read b < adaptive/continue

let i=1
if [ $b == 'ok' ] ; then
	while [ $b == 'ok' ] 
    do
   		
		cd adaptive
			#triangle will automatically add 1 to the file number creating newmesh.(i+1).poly
        	../triangle/./triangle -co2Qq$ang mesh.$i.poly
   		cd ..
			
		let i=$i+1
		#echo 'Removing first line of adaptive/newmesh.'$i'.ele and creating nssim/data_elements.in'
        #Remove the first line of the files and rename
        sed '1d' adaptive/mesh.$i.ele > nssim/elements.in
        sed '1d' adaptive/mesh.$i.node > nssim/nodes.in		
		
		cd adaptive
            if [ $plot == 'T' ] ; then
                ../triangle/./showme mesh.$i.ele
            fi
   		cd ..

		mv adaptive/mesh.$i.poly datafiles
   		mv adaptive/mesh.$i.node datafiles
		mv adaptive/mesh.$i.ele datafiles
		#mv adaptive/mesh.$i.off datafiles
        
		
		cd nssim
		    let j=$i-1
            echo;echo 'Calling Simulation--------------------------------------------------------refinement' $j
            echo $i >> ../summary.out
            ./nssim elements.in nodes.in
            cat summary.out >> ../summary.out
            cat time.out >> ../time.out
            echo
		cd ..

		cd adaptive
            ./adapt ../nssim/nodes.out ../nssim/elements.out \
            ../nssim/solution.out ../nssim/corner.out mesh.$i.poly read.table.o
            cat summary.out >> ../summary.out
		cd ..
		
		read b<adaptive/continue
	done
fi

echo;echo Cleaning up... 
mv adaptive/*.poly datafiles
rm adaptive/summary.out
rm adaptive/refinefile
mv nssim/*.out datafiles/
rm meshgen/domain.out 
rm nssim/simfile
rm nssim/*.in
rm adaptive/continue
rm adaptive/domain
rm adaptive/holefile
rm nssim/domain


echo Generating plots...
cd vectorplot
./vp ../datafiles/nodes.out ../datafiles/velocity.out > plotinfo.out
echo Plot info saved to vectorplot/plotinfo.out
cd ..
mv datafiles/velocity_dir.eps plots
mv datafiles/velocity_vec.eps plots
echo Created plots/velocity_dir.eps and plots/velocity_vec.eps

echo;echo 'Check a summary of the refinement in out summary.out'

open plots/velocity_vec.eps
#geomview datafiles/mesh$i.off


echo;echo

