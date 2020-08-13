// ************************* 
// Tushar Nahar
// File: ShortestPath.cpp
// Version: 12.1.7

// *************************
// All pair Shortest Path Problem: In graph theory, the shortest path problem is the problem of finding a path 
//between two vertices (or nodes) in a graph such that the sum of the weights of its constituent edges is minimized.

// *************************
//MODEL
//This model is run multiple times for diffrent values of source and destination 

// *************************
// Parameters 
//The parameters of the model are
//arc_costs[i][j] 			is the length of the arc from node i and ends at node j

// *************************
//The decision variable 
//arc[i][j] 				is 1 when arc from node i to node j is a part of the shortest path, 0 Otherwise

// *************************
// Objective function
//The objective function of the model is to Minimize the total length of the  path 
//from source node to destination the node.
//Minimize
//Sum(sum over i,j where i not equal to j) 		arc[i][j] x arc_costs [i][j]

//**************************
// Constraints
//The constraints of the model are
//for all intermediate nodes except the source node and destination node for all nodes j

//(sum( sum over i) 			arc[i][j])-(sum( sum over i) arc[j][i])=0
// for source node s		 	sum(sum over j) arc [s][j] =1
//for destination node t		sum(sum over j) arc [j][t] =1

// *************************
// C++ Scripts
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray<IloBoolVarArray>		IloBoolVarArray2;

int main(void) {

	// Open text file to write the results
	ofstream results;
	results.open("output.txt");

	// Setting up cplex environment
	IloEnv env;

	// solution block
	try {
		// Load data
		const char* filename = "input.txt";
		// store the file name in const char filename
		ifstream file(filename);
		if (!file) {
			cerr << "No such file: " << filename << endl;
			throw(-1);
		}
		// Declaring parameters
		IloNumArray2 arc_costs(env);

		// Assign values to parameters
		file >> arc_costs;

		int num_Nodes = arc_costs[0].getSize();

		//Declaring objects
		IloInt i, j, k, s_node, t_node;

		//Declaring variables
		IloBoolVarArray2 arc(env);
		for (i = 0; i < num_Nodes; i++) {
			arc.add(IloBoolVarArray(env, num_Nodes));
		}

		// Checking the input data *********************

		//for (i = 0; i < num_Nodes; i++) {
			//for (j = 0; j < num_Nodes; j++) {
				//env.out() << '\t' << num_Nodes
					//<< '\t' << i
					//<< '\t' << j
					//<< '\t' << arc_costs[i][j] << endl;
			//}
		//}		***************************


		// Solving in for loop

		// Iterating source and sink
		for (s_node = 0; s_node < num_Nodes - 1; s_node++) {
			for (t_node = s_node + 1; t_node < num_Nodes; t_node++) {

				// Create model for particular s and t

				// Declaring objective function
				IloExpr Distance(env);
				for (i = 0; i < num_Nodes; i++) {
					for (j = 1; j < num_Nodes; j++) {
						Distance += arc_costs[i][j] * arc[i][j];
					}
				}
				IloModel mod(env);
				mod.add(IloMinimize(env, Distance));
				Distance.end();


				// Declaring constraints

				// Source flow constraint
				for (i = 0; i < num_Nodes; i++) {
					if (i == s_node) {
						IloExpr s_flow_out(env);
						for (j = 0; j < num_Nodes; j++) {
							s_flow_out += arc[i][j];
						}
						mod.add(s_flow_out == 1);
						s_flow_out.end();
						IloExpr s_flow_in(env);
						for (j = 0; j < num_Nodes; j++) {
							s_flow_in += arc[j][i];
						}
						mod.add(s_flow_in == 0);
						s_flow_in.end();

					}
					// Sink flow constraint
					else if (i == t_node) {
						IloExpr t_flow_in(env);
						for (j = 0; j < num_Nodes; j++) {
							t_flow_in += arc[j][i];
						}
						mod.add(t_flow_in == 1);
						t_flow_in.end();
						IloExpr t_flow_out(env);
						for (j = 0; j < num_Nodes; j++) {
							t_flow_out += arc[i][j];
						}
						mod.add(t_flow_out == 0);
						t_flow_out.end();
					}
					// Flow balance constraints for every node other than source and sink
					else {
						IloExpr flow(env);
						for (j = 0; j < num_Nodes; j++) {
							for (k = 0; k < num_Nodes; k++) {
								flow += arc[j][i] - arc[i][k];
							}
						}
						mod.add(flow == 0);
						flow.end();
					}
				}

				for (i = 0; i < num_Nodes; i++) {
					for (j = 0; j < num_Nodes; j++) {
						if (arc_costs[i][j] == -1) {
							mod.add(arc[i][j] == 0);
						}
					}

				}

				// Solving the problem for particular s & t

				IloCplex cplex(mod);
				cplex.exportModel("output.lp");

				cplex.solve();

				results << '\t' << "Starting Node = " << s_node + 1 << '\t' << "Terminating Node = " << t_node + 1
					<< '\t' << "SPP Distance = " << cplex.getObjValue() << '\t' << "Path = ";


				for (int p = 0; p < num_Nodes; p++) {
					for (int q = 0; q < num_Nodes; q++) {
						if ((cplex.getValue(arc[p][q]) == 1) && (q != t_node)) {
							results << " --> " << q + 1;
						}
					}
				}

				results << endl;
				// ********************
				// Checking output ****
				// ********************
				//env.out() << endl << "SPP Distance = " << cplex.getObjValue() << endl;
				//int p, q;
				//for (p = 0; p < num_Nodes; p++) {
					//for (q = 0; q < num_Nodes; q++) {
						//env.out() << '\t' << cplex.getValue(arc[p][q]) << endl;
					//}
				//}
				// ********************
				mod.end();
			}
		}


	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error: Unknown exception caught!" << endl;
	}

	env.end();

	return 0;
} // END main