add_executable("wave" wave.cc wave.ini wavefem.hh driver.hh bctype.hh)
target_link_dune_default_libraries("wave")

dune_symlink_to_source_files(FILES wave.ini)
