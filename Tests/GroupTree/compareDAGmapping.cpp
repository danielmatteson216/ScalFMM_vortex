
// Keep in private GIT
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <thread>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <omp.h>
using namespace std;

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
enum TaskType {P2M = 0, M2M, M2L, M2L_OUT, L2L, L2P, P2P, P2P_OUT};
string taskNames[] = {"P2M", "M2M", "M2L", "M2L_out", "L2L", "L2P", "P2P", "P2P_out"};
struct Task
{
	TaskType type;
	long long int uniqueId;
	vector<long long int> id;
	int mpiNode;
	int level;
	double perf;
	bool operator==(const Task & other) const
	{
		if(type != other.type || id.size() != other.id.size())
			return false;
		for(size_t i = 0; i < id.size(); ++i)
			if(id[i] != other.id[i])
				return false;
		return true;
	}
	bool operator!=(const Task & other) const
	{
		return !((*this)==other);
	}
	void print(void)
	{
		cout << taskNames[type];
		for(size_t i = 0; i < id.size(); ++i)
			cout << ", " << id[i];
		cout << "(mpi " << mpiNode << ")";
		cout << endl;
	}
};
//Spécialisation de hash pour le type Task
namespace std {
  template <> struct hash<Task>
  {
    size_t operator()(const Task & x) const
    {
		return x.uniqueId;
    }
  };
}
struct DagData
{
	unordered_set<Task> allTask;
	unordered_map<long long int, double> performence;
	int treeHeight;
};

bool parseLine(DagData & dagData, deque<string> & lineElements)
{
	if(lineElements.size() < 1)
		return false;
	Task task;
	if(lineElements.size() >= 14 && lineElements[0] == "P2P_out")
	{
		task.type = P2P_OUT;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(4);
		task.id[0] = stoll(lineElements[9]);
		task.id[1] = stoll(lineElements[10]);
		task.id[2] = stoll(lineElements[11]);
		task.id[3] = stoll(lineElements[12]);
		task.mpiNode = stoi(lineElements[13]);
		task.level = dagData.treeHeight - 1;
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 10 && lineElements[0] == "P2P")
	{
		task.type = P2P;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(4);
		task.id[0] = stoll(lineElements[5]);
		task.id[1] = stoll(lineElements[6]);
		task.id[2] = stoll(lineElements[7]);
		task.id[3] = stoll(lineElements[8]);
		task.mpiNode = stoi(lineElements[9]);
		task.level = dagData.treeHeight - 1;
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 10 && lineElements[0] == "M2L" )
	{
		task.type = M2L;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[5]);
		task.id[2] = stoll(lineElements[6]);
		task.id[3] = stoll(lineElements[7]);
		task.id[4] = stoll(lineElements[8]);
		task.mpiNode = stoi(lineElements[9]);
		task.level = (int)task.id[0];
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 13 && lineElements[0] == "M2L_out")
	{
		task.type = M2L_OUT;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[8]);
		task.id[2] = stoll(lineElements[9]);
		task.id[3] = stoll(lineElements[10]);
		task.id[4] = stoll(lineElements[11]);
		task.mpiNode = stoi(lineElements[12]);
		task.level = (int)task.id[0];
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 13 && lineElements[0] == "M2M")
	{
		task.type = M2M;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[8]);
		task.id[2] = stoll(lineElements[9]);
		task.id[3] = stoll(lineElements[10]);
		task.id[4] = stoll(lineElements[11]);
		task.mpiNode = stoi(lineElements[12]);
		task.level = (int)task.id[0];
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 13 && lineElements[0] == "L2L")
	{
		task.type = L2L;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[8]);
		task.id[2] = stoll(lineElements[9]);
		task.id[3] = stoll(lineElements[10]);
		task.id[4] = stoll(lineElements[11]);
		task.mpiNode = stoi(lineElements[12]);
		task.level = (int)task.id[0];
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 8 && lineElements[0] == "L2P")
	{
		task.type = L2P;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(2);
		task.id[0] = stoll(lineElements[5]);
		task.id[1] = stoll(lineElements[6]);
		task.mpiNode = stoi(lineElements[7]);
		task.level = dagData.treeHeight - 1;
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 8 && lineElements[0] == "P2M")
	{
		task.type = L2P;
		task.uniqueId = stoll(lineElements[1]);
		task.perf = dagData.performence.count(task.uniqueId) == 1 ? dagData.performence[task.uniqueId] : 0.0;
		task.id.resize(2);
		task.id[0] = stoll(lineElements[5]);
		task.id[1] = stoll(lineElements[6]);
		task.mpiNode = stoi(lineElements[7]);
		task.level = 0;
		dagData.allTask.insert(task);
	}
	else
	{
		cout << "No match for " << lineElements[0] << " - " << lineElements.size() << endl;
		return false;
	}
	return true;
}
void split(string& line, string const& delim, deque<string> &result)
{
	size_t prevPos = 0;
	size_t pos = 0;
	while((pos = line.find(delim, prevPos)) != string::npos)
	{
		result.push_back(line.substr(prevPos, pos-prevPos));
		prevPos = pos+delim.size();
	}
	if(prevPos < line.size())
		result.push_back(line.substr(prevPos));
}
bool fillDagData(const char* const filename, DagData & dagData)
{
	std::ifstream fichier(filename, ios::in);  // on ouvre le fichier en lecture

	if(!fichier)  // si l'ouverture a réussi
	{       
		cerr << "Couldn't open " << filename << endl;
		return false;
	}
	string line;
	string delim(", ");
	deque<string> splitLine;
	while(!fichier.eof())
	{
		getline(fichier, line);
		splitLine.clear();
		split(line, delim, splitLine);
		parseLine(dagData, splitLine);
	}
	fichier.close();  // on ferme le fichier
	return true;
}
bool getTaskPerf(string line, double &perf)
{
	perf = stod(line);
	return true;
}
bool getTaskId(string line, long long int & taskId)
{
	size_t pos = line.rfind('_');
	if(pos == string::npos)
		return false;
	taskId = stoll(line.substr(pos));
	return true;
}
void fillPerformanceData(const char* const filename, DagData & dagData)
{
	std::ifstream fichier(filename, ios::in);  // on ouvre le fichier en lecture

	if(!fichier)  // si l'ouverture a réussi
		cerr << "Couldn't open " << filename << endl;

	string line;
	string delim(", ");
	deque<string> splitLine;
	bool getPerf = false;
	long long int taskId;
	double perf;
	while(!fichier.eof())
	{
		getline(fichier, line);
		if(line.size() > 0 && line[0] == 'N')
		{
			if(getTaskId(line.substr(3), taskId))
				getPerf = true;
		}
		else if(getPerf && line.size() > 0 && line[0] == 'S')
		{
			if(getTaskPerf(line.substr(3), perf))
			{
				getPerf = false;;
				dagData.performence[taskId] = perf;
			}
		}
		splitLine.clear();
		split(line, delim, splitLine);
		parseLine(dagData, splitLine);
	}
	fichier.close();  // on ferme le fichier
}
void compareDag(DagData const& dag1, DagData const& dag2, int const treeHeight)
{
	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			long long int notFoundCount[omp_get_num_threads()][treeHeight];
			long long int differenceMapping[omp_get_num_threads()][treeHeight];
			long long int taskCount[omp_get_num_threads()][treeHeight];
			for(int i = 0; i < omp_get_num_threads(); ++i)
			{
				for(int j = 0; j < treeHeight; ++j)
				{
					notFoundCount[i][j] = 0;
					taskCount[i][j] = 0;
					differenceMapping[i][j] = 0;
				}
			}
			for(Task task : dag1.allTask)
			{
				#pragma omp task default(shared) firstprivate(task)
				{
					bool found = false;
					Task sameTask[2];
					int sameTaskId = 0;
					//Count task per level
					if(task.level < treeHeight)
						++taskCount[omp_get_thread_num()][task.level];
					//Look into the second dag to find task in it
					//We look for two task because it may be symetrized
					for(auto it = dag2.allTask.begin(); it != dag2.allTask.end(); ++it)
					{
						if(task == *it)
						{
							sameTask[sameTaskId++] = *it;
							found = true;
							if(sameTaskId == 2)
								break;
						}
					}
					//If the task it not found, print it and count it
					if(found == false)
					{
						#pragma omp critical
						task.print();
						if(task.level < treeHeight)
							++notFoundCount[omp_get_thread_num()][task.level];
					}
					else //else check the mapping
					{
						bool sameNode = false;
						for(int i = 0; i < sameTaskId; ++i)
							if(sameTask[i].mpiNode == task.mpiNode)
								sameNode = true;
						//The tasks are not mapped on the same node	
						if(!sameNode)
							if(task.level < treeHeight)
								++differenceMapping[omp_get_thread_num()][task.level];
					}
				}
			}
			#pragma omp taskwait
			//Sum results by level and print it
			for(int i = 0; i < treeHeight; ++i)
			{
				long long int sum = 0;
				long long int sumDiffMapping = 0;
				long long int sumTaskCount = 0;
				for(int j = 0; j < omp_get_num_threads(); ++j)
					if(taskCount[j][i] > 0)
					{
						sum += notFoundCount[j][i];
						sumDiffMapping += differenceMapping[j][i];
						sumTaskCount += taskCount[j][i];
					}
				//Print only if there is error repported on the level
				if(sum > 0 || sumDiffMapping > 0)
					std::cout << "Diff lvl " << i << " -> " << sum << " (Mapping error : " << sumDiffMapping << "/" << sumTaskCount << ")" << std::endl;
			}
		}
	}
}
int main(int argc, char* argv[])
{
    const FParameterNames ExplicitTrace {
        {"-E"},
        "Simgrid trace from explicit mpi"
    };
	const FParameterNames ImplicitTrace {
		{"-I"} ,
		"Simgrid trace from implicit mpi"
	};
    const FParameterNames Explicit {
        {"-e"},
        "Simgrid trace from explicit mpi"
    };
	const FParameterNames Implicit {
		{"-i"} ,
		"Simgrid trace from implicit mpi"
	};
    const FParameterNames TreeHeight {
        {"-h"},
        "Height of the tree"
    };
    FHelpDescribeAndExit(argc, argv, "Compare DAG mapping", Explicit, Implicit, ExplicitTrace, ImplicitTrace, TreeHeight);

    // Get params
    const char* const explicitFilename = FParameters::getStr(argc,argv,Explicit.options, "scalfmm_explicit.out");
    const char* const implicitFilename = FParameters::getStr(argc,argv,Implicit.options, "scalfmm_implicit.out");
    const int treeHeight = FParameters::getValue(argc,argv,TreeHeight.options, 5);

	DagData implicitData, explicitData;
	bool implicitGood, explicitGood;
	//Read data from files and put it into the DAG data structure
	std::thread explicitThread([&](){
		explicitData.treeHeight = treeHeight;
		explicitGood = fillDagData(explicitFilename, explicitData);
		});
	std::thread implicitThread([&](){
		implicitData.treeHeight = treeHeight;
		implicitGood = fillDagData(implicitFilename, implicitData);
		});
	implicitThread.join();
	explicitThread.join();
	//If every thing went good, start comparing
	if(implicitGood && explicitGood)
	{
		cout << explicitData.allTask.size() << " tasks in explicit." << endl;
		cout << implicitData.allTask.size() << " tasks in implicit." << endl;
		//compareDag(implicitData, explicitData, treeHeight);
		compareDag(explicitData, implicitData, treeHeight);
	}
    return 0;
}
