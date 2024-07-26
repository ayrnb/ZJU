#include "uncertain-core.h"
#include "probability-change-core.h"

std::string getParentDirectory(const std::string& filePath) {
	size_t found = filePath.find_last_of("/\\");
	return filePath.substr(0, found);
}

void load_cores(const char* str, vec_i& cores)
{
	const char* txt = strstr(str, "txt");
	const char* bin = strstr(str, "bin");
	if (txt == NULL && bin == NULL) {
		printf("Usage: file format \"*.txt\" or \"*.bin\"\n");
		exit(1);
	}

	int n = 0;
	FILE* in = NULL;
	if (txt != NULL) in = fopen(str, "r");	//���ı���ȡģʽ���ļ�
	else in = fopen(str, "rb");		//�Զ�����ģʽ���ļ�

	if (in == NULL) {
		printf("No such file: %s\n", str);
		exit(1);
	}
	if (txt != NULL) {
		int x = fscanf(in, "%d", &n);	//x��ʾ�ɹ���ȡ������
		printf("file=%s, n=%d\n", str, n);
		cores.resize(n);
		for (int i = 0; i < n; ++i)
			x = fscanf(in, "%d", &cores[i]);
	}
	else {
		size_t x = fread(&n, sizeof(int), 1, in);
		printf("file=%s, n=%d\n", str, n);
		cores.resize(n);
		for (int i = 0; i < n; ++i)
			x = fread(&cores[i], sizeof(int), 1, in);
	}
	fclose(in);
}

void prinf_core(const char* str, vec_i& cores)	//���ȴ洢�������n��Ȼ�����δ洢c(i)
{
	const char* txt = strstr(str, "txt");
	const char* bin = strstr(str, "bin");
	if (txt == NULL && bin == NULL) {
		printf("Usage: file format \"*.txt\" or \"*.bin\"\n");
		exit(1);
	}
	int n = cores.size();
	FILE* in = NULL;
	if (bin != NULL) in = fopen(str, "wb");	//�Զ�����д��ģʽ���ļ�
	else in = fopen(str, "w");		//���ı�д��ģʽ���ļ�

	if (in == NULL) {
		printf("No such file: %s\n", str);
		exit(1);
	}
	if (bin != NULL) {
		fwrite(&n, sizeof(int), 1, in);
		for (int i = 0; i < n; ++i)
			fwrite(&cores[i], sizeof(int), 1, in);
	}
	else {
		//printf("n=%d\n", n);
		fprintf(in, "%d\n", n);		//���ı���ʽд���ļ�
		for (int i = 0; i < n; ++i)
			fprintf(in, "%d\n", cores[i]);
	}
	fclose(in);
}

//��ͼ�洢��bin�ļ���--��ȡ�ٶȸ���
void printf_thres(const std::string& filePath, const std::vector<std::vector<double> >& thres) {
	size_t txtPos = filePath.find(".txt");
	size_t binPos = filePath.find(".bin");
	if (txtPos == std::string::npos && binPos == std::string::npos) {
		printf("Usage: file format \"*.txt\" or \"*.bin\"\n");
		exit(1);
	}
	int k = thres.size();
	int n = thres[0].size();
	FILE* in = NULL;
	if (binPos != std::string::npos) in = fopen(filePath.c_str(), "wb");    //�Զ�����д��ģʽ���ļ�
	else in = fopen(filePath.c_str(), "w");        //���ı�д��ģʽ���ļ�

	if (in == NULL) {
		printf("No such file: %s\n", filePath.c_str());
		exit(1);
	}
	if (binPos != std::string::npos) {
		fwrite(&k, sizeof(int), 1, in);
		//fwrite(&n, sizeof(int), 1, in);
		for (int i = 0; i < k; i++) {
			//fwrite(&i, sizeof(int), 1, in);
			for (int j = 0; j < n; j++) {
				fwrite(&thres[i][j], sizeof(double), 1, in);
			}
		}
	}
	else {
		for (int i = 0; i < k; i++) {	//txt�ļ� - �ɶ��Ը�ǿ
			fprintf(in, "%d :\n", i);
			for (int j = 0; j < n; j++) {
				fprintf(in, "%d :", j);
				fprintf(in, " %lf\n", thres[i][j]);
			}
		}
	}
	fclose(in);
}


//ʹ��bin�ļ���ȡthres
void load_thres(const std::string& filePath, std::vector<std::vector<double> >& thres) {
	double tm = omp_get_wtime();
	size_t binPos = filePath.find(".bin");
	if (binPos == std::string::npos) {
		printf("Usage: file format \"*.bin\"\n");
		exit(1);
	}

	int kk = 0;
	FILE* in = NULL;
	in = fopen(filePath.c_str(), "rb");
	if (in == NULL) {
		printf("No such file: %s\n", filePath.c_str());
		exit(1);
	}

	size_t x = fread(&kk, sizeof(int), 1, in);
	cout << "read kmax:" << kk << endl;
	for (int i = 0; i < thres.size(); i++) {
		for (int j = 0; j < thres[0].size(); j++) {
			x = fread(&thres[i][j], sizeof(double), 1, in);
		}
	}
	tm = omp_get_wtime() - tm;
	std::cout << "��ȡ��ʼthres�ļ�---- time��" << tm << endl;
}



//����ߵĳ�ʼ�����㷨
vector<vector<double> > insertEdges(string infile, string parentPath, double scale) {
	Uncertain_Core uc;

	//��ԭ�����ļ� ���ѡ�� ��ʼ�� ������ʼͼ
	uc.edge_selected_bin(infile, scale);	//scale - unselected

	int n = uc.get_nm();
	uc.get_core();
	int kmax = uc.get_kmax();
	cout << "kmax:" << kmax << endl;
	vector<vector<double> > thres(kmax, vector<double>(n));
	uc.Initial_threshold_compute_map(thres);

	//std::string thresInFile = parentPath + "/datas/decomposition/d-initial-Flickr.bin";	//��ȡ�ֽ��ļ�
	//load_thres(thresInFile, thres);	

	//ȷ��ͼ�����
	//uc.insert_core_compare();
	uc.insert_threshold_compare(thres);
	return thres;
}


//ɾ���ߵĳ�ʼ�����㷨
void deleteEdges(string readfile, string parentPath, double scale) {
	Uncertain_Core uc;
	uc.read_bin(readfile);
	int n = uc.get_nm();
	uc.get_core();
	int kmax = uc.get_kmax();
	vector<vector<double> > thres(kmax, vector<double>(n));	//��ʼ���� ���� ��洢��then��ȡ����

	//��ʼ����thres	
	uc.Initial_threshold_compute_map(thres);
	//std::string thresInFile = parentPath + "/maintence datas/d-initial-try.bin";
	//printf_thres(thresInFile, thres);

	//��ȡthres	
	//load_thres(thresInFile, thres);

	uc.delete_threshold_compare(thres, scale);
	//uc.delete_compare_range(thres, scale);

	/*std::string outfile = parentPath + "/maintence datas/d-read-try.txt";
	printf_thres(outfile, thres);*/
}


int main()
{
	// ��ȡ main.cpp ��·����main����Ŀ¼�������ļ������main.cpp �����·�������������ļ��ľ���·��
	std::string currentPath = __FILE__;
	std::string parentPath = getParentDirectory(currentPath);
	std::string dataFilePath = "datas/Fruit-Fly.bin";	
	std::string infile = parentPath + "/" + dataFilePath;

	double scale = 0.8;		//ѡ�������
	insertEdges(infile, parentPath, scale);	
	deleteEdges(infile, parentPath, scale);
	return 0;
	
	
	Probability_Core pc;
	//pc.creat_bin(infile);	//����bin�ļ�

	pc.read_bin(infile);
	int n = pc.get_nm();

	double tm = omp_get_wtime();
	pc.get_core();
	int kmax = pc.get_kmax();
	vector<vector<double> > thres(kmax, vector<double>(n));	//��ʼ���� ���� ��洢��then��ȡ����
	pc.Initial_threshold_compute_map(thres);
	tm = omp_get_wtime() - tm;
	std::cout << "Time:" << tm << endl;


	pc.increase_threshold_compare(thres, 0.8, 1000);	//��Ե��������
	pc.decrease_threshold_compare(thres, 0.8, 1000);	//��Ե���ʼ�С


	//save decomposition result
	/*std::string outfile = parentPath + "/datas/d-Flickr.txt";
	printf_thres(outfile, thres);*/


	system("pause");
	return 0;
}