add_library_unity(gpgomea_genotype OBJECT Node.cpp SingleNode.cpp Multitree.cpp)
set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:gpgomea_genotype>
    PARENT_SCOPE)

