file(REMOVE_RECURSE
  "vmc_parallel.pdb"
  "vmc_parallel"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/vmc_parallel.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
