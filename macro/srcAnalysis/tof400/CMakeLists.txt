set(INCLUDE_DIRECTORIES
		${BASE_INCLUDE_DIRECTORIES}
		${CMAKE_SOURCE_DIR}/bmndata
		${CMAKE_SOURCE_DIR}/bmndata/cbm
		${CMAKE_SOURCE_DIR}/bmndata/math
		${CMAKE_SOURCE_DIR}/bmnfield
		${CMAKE_SOURCE_DIR}/silicon
		${CMAKE_SOURCE_DIR}/alignment
		${CMAKE_SOURCE_DIR}/globaltracking
		${CMAKE_SOURCE_DIR}/uni_db
		${CMAKE_SOURCE_DIR}/uni_db/db_classes
		${CMAKE_SOURCE_DIR}/steering
		${CMAKE_SOURCE_DIR}/gem
		${CMAKE_SOURCE_DIR}/tof1
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

link_directories(${LINK_DIRECTORIES})

	set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin/tofCalib")

	add_executable(gammaPeak gammaPeak.cpp)
target_link_libraries(gammaPeak BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

	add_executable(multiplicity multiplicity.cpp)
target_link_libraries(multiplicity BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

	add_executable(alignStrips stripAlignment/alignStrips.cpp)
target_link_libraries(alignStrips BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

	add_executable(tofTimeWalk timeWalk/tofTimeWalk.cpp)
target_link_libraries(tofTimeWalk BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)

	add_executable(tofCluster clustering/tofCluster.cpp)
target_link_libraries(tofCluster BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
	
	add_executable(fastTimeZ2 fastTimeZ2.cpp)
target_link_libraries(fastTimeZ2 BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
	
	add_executable(tofAnalysis tofAnalysis.cpp)
target_link_libraries(tofAnalysis BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
	
	add_executable(class_tofAnalysis class_tofAnalysis.cpp)
target_link_libraries(class_tofAnalysis BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
	
	add_executable(testUnion testUnion.cpp)
target_link_libraries(testUnion BmnData Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
