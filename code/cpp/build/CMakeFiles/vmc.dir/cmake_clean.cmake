file(REMOVE_RECURSE
  "vmc.pdb"
  "vmc"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/vmc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
