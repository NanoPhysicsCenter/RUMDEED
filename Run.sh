!#bash.rc
make clean all
if [ -f build/Vacuum-MD.out ]
then
    cp build/Vacuum-MD.out data
    cd data
    time ./Vacuum-MD.out
else
    echo 'Compilation Failed'
fi
