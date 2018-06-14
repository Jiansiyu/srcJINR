set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/bmndata
 ${CMAKE_SOURCE_DIR}/globaltracking
 ${CMAKE_SOURCE_DIR}/gem
 ${CMAKE_SOURCE_DIR}/silicon
 ${CMAKE_SOURCE_DIR}/uni_db
 ${CMAKE_SOURCE_DIR}/uni_db/db_classes
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${Boost_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIR}
  )

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
 ${GSL_LIBRARY_DIRS}
)


set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin/timeWalk")

link_directories(${LINK_DIRECTORIES})

add_executable(skimForNonT0 skimForNonT0.cpp)
target_link_libraries(skimForNonT0 BmnData Silicon TOF1 UniDb Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

add_executable(skimForT0 skimForT0.cpp)
target_link_libraries(skimForT0 BmnData Silicon TOF1 UniDb Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

add_executable(nonT0-timewalk nonT0-timewalk.cpp)
target_link_libraries(nonT0-timewalk BmnData Silicon TOF1 UniDb Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

add_executable(T0_timewalk T0_timewalk.cpp)
target_link_libraries(T0_timewalk BmnData Silicon TOF1 UniDb Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
