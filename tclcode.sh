#! /bin/sh

echo "/* Automatically generated. Do not edit ! */"
echo "#include \"tclutil.h\""

for i in $*; do
  j=`echo $i | sed -e 's/\.tcl//'`
  sh tcl2c.sh ${j}_tcl $i
  lst="$lst
{\"$j\", ${j}_tcl},"
done

echo "TCLCODE tclcode_pointers[] = {
$lst
{\"\", (char*)0}
};"

