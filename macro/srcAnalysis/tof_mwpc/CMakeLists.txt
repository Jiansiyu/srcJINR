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
		${CMAKE_SOURCE_DIR}/mwpc
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

	set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin/fullEvent")

	
	add_executable(tof_mwpc tof_mwpc.cpp)
target_link_libraries(tof_mwpc BmnData BmnField Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
	
	add_executable(playSegments playSegments.cpp)
target_link_libraries(playSegments BmnData BmnField Core Base Sts Silicon TOF1 UniDb BmnSteering Gem GlobalTracking KF Mwpc Dch ${ROOT_LIBRARIES} XMLParser Minuit Spectrum GenVector MathMore)
