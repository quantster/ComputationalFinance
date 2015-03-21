INC_DIR:= ./include
SRC_DIR:= ./source
SRC:= $(SRC_DIR)/BlackScholesFormulas.cpp $(SRC_DIR)/Normals.cpp $(SRC_DIR)/Random1.cpp $(SRC_DIR)/Vanilla1.cpp $(SRC_DIR)/PayOff2.cpp MonteCarloSimulator.cpp 
hw1 : $(SRC) HW1main.cpp
	g++ HW1main.cpp $(SRC) -I$(INC_DIR) -o hw1 -wall
