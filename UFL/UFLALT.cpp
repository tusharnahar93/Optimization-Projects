
#include<ilcplex/ilocplex.h>
ILOSTLBEGIN
typedef IloArray<IloNumVarArray> IloNumVarArray2;
int main(int argc, char **argv)
{
	ofstream myfile;
	myfile.open("output.txt");
	IloEnv env;
	try {
		const char* filename = "UFL.dat";
		if (argc >= 2) filename = argv[1];
		ifstream file(filename);
		if (!file) {
			cerr << "No_such_file:_" << filename << endl;
			throw(-1);
		}

		//PARAMETERS 
		IloNumArray demand(env);  /*Declaring parameter called demand for the client demand*/
		IloNumArray unitcost(env); /*Declaring parameter called unitcost for the cost per unit*/
		IloNumArray fixedcost(env); /*Declaring parameter called fixed cost for the fixed cost of building a facility*/

		//IloNumArray2 adjnew(env);
		file >> demand >> unitcost >> fixedcost;  /* Reads the data from file and stores the values in order */
		IloInt nq = demand.getSize(); /*Parameter q stores the number of clients*/
		IloInt np = fixedcost.getSize(); /*Parameter p stores the number of facilities*/



		// VARIABLES

		IloModel mod(env); /*Creating two separate model variables for the original MIP and for the LP Relaxation*/
		IloModel relax(env);

		IloNumVarArray2 X(env); //creating a two dimensional variable X which denotes the portion of the demand of client q satsfied by facility p.
		for (int p = 0; p < np; p++)
		{
			X.add(IloNumVarArray(env, nq, 0, IloInfinity));
		}

		IloNumVarArray Y(env); //creating a two dimensional variable Y in the environment.
		for (int p = 0; p < np; p++)
		{
			Y.add(IloNumVar(env, 0, 1, ILOBOOL)); // Y is a binary variable. Y=1 if the production happens at facility p, else 0.
		}
		//OBJECTIVE is to minimize the total cost in the path
		IloExpr TotalUnitCost(env); //Creating Expressions for the Costs
		IloExpr TotalFixedCost(env);
		IloExpr TotalCost(env);
		for (int p = 0; p < np; p++)
		{
			TotalFixedCost += fixedcost[p] * Y[p]; //Calculation of fixed costs
			for (int q = 0; q < nq; q++)
				TotalUnitCost += unitcost[(p)* nq + q] * X[p][q]; //Calculation of Unit Distribution costs
		}
		TotalCost = TotalFixedCost + TotalUnitCost; // The total cost expression is sum of the Fixed Costs and Unit distribution costs 
		mod.add(IloMinimize(env, TotalCost)); //Defining the objective for the MIP Problem
		relax.add(IloMinimize(env, TotalCost)); //Defining the objective for the LP Problem
		relax.add(IloConversion(env, Y, ILOFLOAT));// This line converts the binary variable into a float type (Relaxes)
		TotalFixedCost.end(); //Ending the created expressions for Costs
		TotalUnitCost.end();
		TotalCost.end();

		//Constraints

		/* */
		for (int q = 0; q < nq; q++)
		{

			IloExpr DemandSat(env); //Creating  an expression for demand satisfied for every client
			for (int p = 0; p < np; p++)
			{
				DemandSat += X[p][q]; // DemandSat is equal to the sum of the demand portions satisfied by each facility for that client
			}
			mod.add(DemandSat == demand[q]); // DemandSat for each client must be eqaul to the demand of that client. This constraint is added to the MIP model.
			relax.add(DemandSat == demand[q]);// The same constraint is added for the LP relaxation model.
			DemandSat.end(); //Ending the expression 

		}

		IloInt TotalDemand = 0; // Creating a integer to store the total demand. Initially assigning it to 0.
		for (int q = 0; q < nq; q++)
		{
			TotalDemand += demand[q]; //Summing the demands of all the clients and storing it in total demand. 
		}
		for (int p = 0; p < np; p++)
		{
			

			for (int q = 0; q < nq; q++)
			{
				
				mod.add(X[p][q]<=demand[q]*Y[p]);/* The demand portion of client q satisfied by facility p must be less than  
				                                    the total demand of client q only if the facility p is open. This constraint is added for every p and q. 
													This constraint is added to the MIP Problem*/

				relax.add(X[p][q] <= demand[q] * Y[p]);// The same constraint is added to the LP relaxation model.

			}
			
			
		}


		IloCplex cplex1(relax); //Creating a cplex instance to solve the relax model
		IloCplex cplex(mod); // Creating a cplex instance to solve the MIP model
		cplex.exportModel("outputAlt.lp");
		cplex1.exportModel("output1Alt.lp");

		cplex1.solve();  //Solving the LP Relaxation
		IloNum LPtime;  // Declating a variable to store the solving time for LP Relaxation

		LPtime = cplex1.getCplexTime(); // LP relaxation solve time is stored in LPtime.

		cplex.setParam(IloCplex::CutPass, -1); //Turning off Default Cuts
		cplex.solve(); // Solve the MIP Problem 

		//////////////////////////////////////          SECTION A         ///////////////////////////////


		env.out() << endl << "List of Production Facilities with Non-Zero Production and their total Production" << endl;

		IloInt sum = 0; //Creating an integer to store the total production
		for (int p = 0; p < np; p++)
		{
			if (cplex.getValue(Y[p]) == 1) //Checking if some production occurs at facility p
			{

				for (int q = 0; q < nq; q++)
				{
					sum += cplex.getValue(X[p][q]); //Adding all the production that occurs at facility p

				}
				env.out() << endl << "Facility" << '\t' << p << '\t' << "TotalProduction" << '\t' << sum << '\t';
				sum = 0;
			}
		}

		env.out() << endl;
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t'<< cplex.getObjValue() << endl;

		env.out() << endl << "Number of Integer variables" << '\t' << np << endl;
		env.out() << endl << "Number of Continuous Variables" << '\t' << np * nq << endl;
		env.out() << endl << "Number of Constraints (Including binary and non negativity constraint)" << '\t' << np + nq + 2* np * nq << endl;

		env.out() << endl << "Run time to solve the LP relaxation Problem " << '\t' << LPtime << '\t' << "seconds" << endl;

		env.out() << endl << "Optimal Objective Value of LP Relaxation" << '\t' << cplex1.getObjValue() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplex.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplex.getNnodes() << endl;

		env.out() << endl << "Percentage Gap" << '\t' << cplex.getMIPRelativeGap() << endl;

		env.out() << endl << "Number of Mixed Integer Rounding Cuts" << '\t' << cplex.getNcuts(cplex.CutMir) << endl;

		env.out() << endl << "Number of Flow Cuts" << '\t' << cplex.getNcuts(cplex.CutFlowCover) << endl;

		env.out() << endl << "Number of Implied Bound cuts" << '\t' << cplex.getNcuts(cplex.CutImplBd) << endl;

		env.out() << endl << "Number of Cuts by CPLEX (sum of the thre cutsmentioned above)" << cplex.getNcuts(cplex.CutMir) + cplex.getNcuts(cplex.CutFlowCover) + cplex.getNcuts(cplex.CutImplBd) << endl;


		


	}

	catch (IloException& e) {
		cerr << "Error:_" << e << endl;
	}
	catch (...) {
		cerr << "Error:_Unknown_exception_caught!" << endl;
	}
	int BP1;
	cin >> BP1;
	env.end();
	return 0;
}