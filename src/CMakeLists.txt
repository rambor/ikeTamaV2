MESSAGE("** CURRENT SOURCE INCLUDE DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND SOURCE_FILES
        "${CMAKE_CURRENT_LIST_DIR}/PointSetModel.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SubUnit.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/KDE.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Component.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Anneal.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Objective.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/DensityMapper.cpp"
        )

MESSAGE("** SOURCE FILES => ${SOURCE_FILES}")