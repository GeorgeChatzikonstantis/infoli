# This useful tool purges the repo from unnecessary files
# that are not tracked by git.

cd ../src
make clean
cd ../tools
make clean
cd ../run
echo ""
echo "-- Cleaning I/O --"
echo ""
cd input
ls | grep -v readme | xargs rm -f
cd ../output
ls | grep -v readme | xargs rm -f
cd ../../scripts
