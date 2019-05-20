#!/bin/bash
## create html with links for output_dir 
## usage:  make_html.sh dataset_path dataset_extra_file_path 
defuse_out=$1
extra_files_path=$2
if [ -e $defuse_out ]
then
  echo '<html><head><title>Defuse Output</title></head><body>' > $defuse_out
  echo '<h2>Defuse Output Files</h2><ul>' >>  $defuse_out
  pushd $extra_files_path
  for f in `find -L . -maxdepth 1 -type f`;
   do fn=`basename ${f}`; echo '<li><a href="'${fn}'">'${fn}'</a></li>' >>  $defuse_out;
  done
  popd
  echo '</ul>' >> $defuse_out
  echo '</body></html>' >>  $defuse_out
fi

