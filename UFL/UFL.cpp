///Tushar Nahar 
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
			IloExpr DemandSat1(env); //Creating an expression DemandSat1 to store the total production in each facility p
			for (int q = 0; q < nq; q++)
			{
				DemandSat1 += X[p][q];// Assigning the total production of facility p to DemandSat1

			}
			mod.add(DemandSat1 <= TotalDemand * Y[p]); /* The total Production in each facility p should be less than the total demand 
													   if the production happens in facility p. This constraint is added to the MIP*/
			relax.add(DemandSat1 <= TotalDemand * Y[p]);// The same constraint is added to the LP Relaxation.

			DemandSat1.end();
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
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t' << cplex.getObjValue() << endl;

		env.out() << endl << "Number of Integer variables" << '\t' << np << endl;
		env.out() << endl << "Number of Continuous Variables" << '\t' << np * nq << endl;
		env.out() << endl << "Number of Constraints (Including binary and non negativity constraint)" << '\t' << np + nq + 2 * np * nq << endl;

		env.out() << endl << "Run time to solve the LP relaxation Problem " << '\t' << LPtime << '\t' << "seconds" << endl;

		env.out() << endl << "Optimal Objective Value of LP Relaxation" << '\t' << cplex1.getObjValue() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplex.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplex.getNnodes() << endl;

		env.out() << endl << "Percentage Gap" << '\t' << cplex.getMIPRelativeGap() << endl;

		env.out() << endl << "Number of Mixed Integer Rounding Cuts" << '\t' << cplex.getNcuts(cplex.CutMir) << endl;

		env.out() << endl << "Number of Flow Cuts" << '\t' << cplex.getNcuts(cplex.CutFlowCover) << endl;

		env.out() << endl << "Number of Implied Bound cuts" << '\t' << cplex.getNcuts(cplex.CutImplBd) << endl;

		env.out() << endl << "Number of Cuts by CPLEX (sum of the thre cutsmentioned above)" << cplex.getNcuts(cplex.CutMir) + cplex.getNcuts(cplex.CutFlowCover) + cplex.getNcuts(cplex.CutImplBd) << endl;

		

		////////////////////////////////////////////                  SECTION B                  /////////////////////////////////////

		IloCplex cplexb(mod); //Creating a cplex instance
		cplexb.exportModel("outputB.lp");
		env.out() << endl << "B) Default Cuts Turned off" << endl;
		cplexb.setParam(IloCplex::CutPass, -1); //Turning off default cuts

		cplexb.solve();

		env.out() << endl << "List of Production Facilities with Non-Zero Production and their total Production" << endl;

		IloInt sum1 = 0; //Creating an integer to store the total production
		env.out() << endl << "List of Production Facilities with Non-Zero Production and their total Production" << endl;

		
		for (int p = 0; p < np; p++)
		{
			if (cplexb.getValue(Y[p]) == 1) //Checking if some production occurs at facility p
			{

				for (int q = 0; q < nq; q++)
				{
					sum1 += cplexb.getValue(X[p][q]); //Adding all the production that occurs at facility p

				}
				env.out() << endl << "Facility" << '\t' << p << '\t' << "TotalProduction" << '\t' << sum1 << '\t';
				sum1 = 0;
			}
		}

		env.out() << endl;
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t' << cplexb.getObjValue() << endl;

		env.out() << endl << "Number of Integer variables" << '\t' << np << endl;
		env.out() << endl << "Number of Continuous Variables" << '\t' << np * nq << endl;
		env.out() << endl << "Number of Constraints (Including binary and non negativity constraint)" << '\t' << np + nq + np + np * nq << endl;

		env.out() << endl << "Run time to solve the LP relaxation Problem " << '\t' << LPtime << '\t' << "seconds" << endl;

		env.out() << endl << "Optimal Objective Value of LP Relaxation" << '\t' << cplex1.getObjValue() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplexb.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplexb.getNnodes() << endl;

		env.out() << endl << "Percentage Gap" << '\t' << cplexb.getMIPRelativeGap() << endl;

		env.out() << endl << "Number of Mixed Integer Rounding Cuts" << '\t' << cplexb.getNcuts(cplexb.CutMir) << endl;

		env.out() << endl << "Number of Flow Cuts" << '\t' << cplexb.getNcuts(cplexb.CutFlowCover) << endl;

		env.out() << endl << "Number of Implied Bound cuts" << '\t' << cplexb.getNcuts(cplexb.CutImplBd) << endl;

		env.out() << endl << "Number of Cuts by CPLEX (sum of the thre cuts mentioned above)" << cplexb.getNcuts(cplexb.CutMir) + cplexb.getNcuts(cplexb.CutFlowCover) + cplexb.getNcuts(cplexb.CutImplBd) << endl;

		/////////////////////////////////////////////               SECTION D, part a)                     /////////////////////////////////////////
		
		
		IloCplex cplexd1(mod); //Creating a cplex instance
		cplexd1.exportModel("outputd1.lp");
		env.out() << endl << "Case a: Preprocessing is turned off" << endl;
		cplexd1.setParam(IloCplex::CutPass, -1); //Default cuts turned off
		cplexd1.setParam(IloCplex::PreInd, 0); //Preprocessing turned off

		cplexd1.solve(); //solving the model with the above settngs

		env.out() << endl << "List of Production Facilities with Non-Zero Production and their total Production" << endl;

		IloInt sum6 = 0; // Creating an integer to store the total production
		for (int p = 0; p < np; p++)
		{
			if (cplexd1.getValue(Y[p]) == 1) //Checking if some production occurs at facility p
			{

				for (int q = 0; q < nq; q++)
				{
					sum6 += cplexd1.getValue(X[p][q]);  //Adding all the production that occurs at facility p

				}
				env.out() << endl << "Facility" << '\t' << p << '\t' << "TotalProduction" << '\t' << sum6 << '\t';
				sum6 = 0;
			}
		}

		env.out() << endl;
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t' << cplexd1.getObjValue() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplexd1.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplexd1.getNnodes() << endl;




		/////////////////////////////////                   SECTION D, part b) Best-Bound Search              //////////////////////////
		env.out() << endl << "Case b: Using MIP Node Selection Strategies" << endl;
		IloCplex cplexd2(mod); //Creating a cplex instance
		cplexd2.exportModel("outputd2.lp");
		env.out() << endl << "Best-bound search" << endl;
		cplexd2.setParam(IloCplex::CutPass, -1); //Turning Off Default Cuts
		cplexd2.setParam(IloCplex::NodeSel, 1);  //Using Best Bound Search
		cplexd2.solve();

		env.out() << endl << "List of Production Facilities with Non-Zero Production and their total Production" << endl;

		IloInt sum2 = 0;
		for (int p = 0; p < np; p++)
		{
			if (cplexd1.getValue(Y[p]) == 1) //Checking if some production occurs at facility p
			{

				for (int q = 0; q < nq; q++)
				{
					sum6 += cplexd1.getValue(X[p][q]);  //Adding all the production that occurs at facility p

				}
				env.out() << endl << "Facility" << '\t' << p << '\t' << "TotalProduction" << '\t' << sum2 << '\t';
				sum2 = 0;
			}
		}

		env.out() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplexd2.getCplexTime() << '\t' << "seconds" << endl;

		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplexd2.getNnodes() << endl;



		/////////////////////////////////             SECTION D, Part b), Depth-First            ////////////////////////////////// 


		env.out() << endl << "Depth - First search" << endl;
		IloCplex cplexd3(mod); //Creating a cplex instance
		cplexd3.exportModel("outputd3.lp");
		cplexd3.setParam(IloCplex::CutPass, -1); // Default cuts turned off
		cplexd3.setParam(IloCplex::NodeSel, 0); //Using Depth First Strategy
		cplexd3.solve();

		env.out() << endl;
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t' << cplexd3.getObjValue() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplexd3.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplexd3.getNnodes() << endl;

		

	  //////////////////////////////////////        SECTION D, Part c),  Maximum Infeasibility         ///////////////////////////////

	    IloCplex cplexd4(mod); //Creating a cplex instance
		cplexd4.exportModel("outputd4.lp");

		env.out() << endl << "Using MIP Variable Selection Strategies" << endl;
		env.out() << endl << "Maximum infeasibility rule" << endl;
		cplexd4.setParam(IloCplex::CutPass, -1); //Default cuts turned off
		cplexd4.setParam(IloCplex::VarSel, 1); //Using maximum infeasibility rule
		cplexd4.solve();
		

		env.out() << endl;
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t' << cplexd4.getObjValue() << endl;

		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplexd4.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplexd4.getNnodes() << endl;


		//////////////////////////////////////        SECTION D, Part c),  Minimum Infeasibility         ///////////////////////////////


		env.out() << endl << "Minimum infeasibility rule" << endl;
		IloCplex cplexd5(mod); //Creating a cplex instance
		cplexd5.exportModel("outputd5.lp");
		cplexd5.setParam(IloCplex::CutPass, -1); //Turning off default cuts
			cplexd5.setParam(IloCplex::VarSel, -1); //Using minimum infeasibilty 
			cplexd5.solve();
	

		env.out() << endl;
		env.out() << endl << "Optimal Objective Value of MIP Instance" << '\t'<< cplexd5.getObjValue() << endl;


		env.out() << endl << "Run time to solve the MIP Problem " << '\t' << cplexd5.getCplexTime() << '\t' << "seconds" << endl;
		env.out() << endl << "Number of branch and bound nodes" << '\t' << cplexd5.getNnodes() << endl;

		



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