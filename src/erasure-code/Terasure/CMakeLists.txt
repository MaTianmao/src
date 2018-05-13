# Terasure plugin

set(Terasure_utils_src
  ErasureCodePluginTerasure.cc
  ErasureCodeTerasure.cc)
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -mavx2")
add_library(Terasure_utils OBJECT ${Terasure_utils_src})
add_dependencies(Terasure_utils ${CMAKE_SOURCE_DIR}/src/ceph_ver.h)


set(Terasure_srcs
  Terasure/src/reed_solomon.c
  Terasure/src/hh.c
  Terasure/src/mbr.c
  Terasure/src/msr.c)
add_library(Terasure_objs OBJECT ${Terasure_srcs}) 

set(ec_Terasure_objs
  $<TARGET_OBJECTS:Terasure_objs>
  $<TARGET_OBJECTS:Terasure_utils>
  $<TARGET_OBJECTS:erasure_code_objs>)

add_library(ec_Terasure SHARED ${ec_Terasure_objs})
set_target_properties(ec_Terasure PROPERTIES
  INSTALL_RPATH "")
target_link_libraries(ec_Terasure ${EXTRALIBS})
install(TARGETS ec_Terasure DESTINATION ${erasure_plugin_dir})

# legacy libraries
foreach(flavor ${Terasure_legacy_flavors})
  set(plugin_name "ec_Terasure_${flavor}")
  add_library(${plugin_name} SHARED ${ec_Terasure_objs})
  set_target_properties(${plugin_name} PROPERTIES
    INSTALL_RPATH "")
  install(TARGETS ${plugin_name} DESTINATION ${erasure_plugin_dir})
  add_dependencies(ec_Terasure ${plugin_name})
endforeach()

if(WITH_EMBEDDED)
  add_library(cephd_ec_Terasure STATIC
    $<TARGET_OBJECTS:Terasure_objs>
    ${Terasure_utils_src})
  set_target_properties(cephd_ec_Terasure PROPERTIES COMPILE_DEFINITIONS BUILDING_FOR_EMBEDDED)
endif()
