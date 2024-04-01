#include "uncertain-core.h"

//��m�����б�Ǵ�ɾ���ı�
void Uncertain_Core::selectRandomEdges(double scale) {
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);
    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //ȷ��ÿ����ֻ����һ��
                if (rand() > scale_size) {  //��ɾ��
                    std::cout << "unselected-vertex:" << i << " nei:" << u << endl;
                    std::pair<int, int> p1 = std::make_pair(i, u);
                    double p = adj[i][j].p;
                    unselected.push_back(make_pair(p1, p));
                }
            }
        }
    }
    std::cout << "unselected edges:" << unselected.size() << endl;
}


//����ɾ���ߣ��Ƚ�coreά���㷨 & ������(����insert���㺯��)
void Uncertain_Core::delete_threshold_compare(vector<vector<double> >& thres, double scale) {
    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;    
    long long int candi_sum = 0;    //ÿ������Ѱ�Һ�ѡ����������ۼƺ�

    selectRandomEdges(scale);

    int um = unselected.size();
    for (int i = 0; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        std::cout << "��i���ߣ�" << i << " u:" << u << " v:" << v << endl;

        //ά��ͼ�ṹ
        bool flag = false;
        int d = deg[u];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[u][j].u;
                if (w == v) {
                    flag = true;
                }
            }
            else
            {
                adj[u][j - 1].u = adj[u][j].u;
                adj[u][j - 1].p = adj[u][j].p;
            }
        }

        flag = false;
        d = deg[v];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[v][j].u;
                if (w == u) {
                    flag = true;
                }
            }
            else
            {
                adj[v][j - 1].u = adj[v][j].u;
                adj[v][j - 1].p = adj[v][j].p;
            }
        }
        deg[u]--;
        deg[v]--;

        //int mincore = min(core[u], core[v]);
        vector<int> candi_core(n);   
        vector<int> recompute_core(n);  //�Ƚ�core�Ƿ���ͬ
        vector<vector<double> > range = delete_threshold_compute_candidate(thres, u, v, tm_1, candi_count);
        //vector<vector<double> > range = delete_threshold_compute_update_range(thres, u, v, tm_1, candi_count);
        //vector<vector<double> > range = delete_threshold_compute_shrink_low(thres, u, v, tm_1, candi_count);   
        //vector<vector<double> > range = delete_threshold_compute_batchUp(thres, u, v, tm_1, candi_count);
        std::cout << "delete_compute_candidate:" << tm_1 << endl;
        candidate_tm += tm_1;
        candi_sum += candi_count;
        for (int t = 0; t < n; t++) {
            candi_core[t] = core[t];
        }

        delete[] core;
		vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);
        std::cout << "delete_threshold_recompute:" << tm_2 << endl;
        recompute_tm += tm_2;
        for (int t = 0; t < n; t++) {
            recompute_core[t] = core[t];
        }

        //�Ƚ�ɾ���ߺ�ı仯��ֵ
        //compareArrays(copy_thres, thres_2, range[0], range[1], mincore - 1);    //���ܱȽ�core��С�ĵ�
        bool compare1 = compareArraysTrueOrFalse(thres, thres_2);
        if (!compare1) {
            std::cout << "�����Ƿ���ͬ-----------------------------------------------------------------------��" << compare1 << endl;
        }
        bool compare2 = compareCoreArrays(candi_core, recompute_core);
        if (!compare2) {
            std::cout << "core�Ƿ���ͬ-----------------------------------------------------------------------��" << compare2 << endl;
        }

        //��newArray���Ƹ�thres
        if (thres.size() == thres_2.size()) {
            thres = thres_2;
        }
        else
        {
            thres.resize(thres_2.size());
            for (size_t i = 0; i < kmax; i++) {
                thres[i] = thres_2[i];
            }
        }        
    }
    //ʱ���Ƚ�+ɾ��unselected
    std::cout << "delete_compute_candidate_time:" << candidate_tm << endl;
    std::cout << "delete_threshold_recompute_time:" << recompute_tm << endl;
    std::cout << "candi sum:" << candi_sum << endl;
    unselected.clear();
}


//��δ����nei��kporb�������� -- �ú��������� δ��ʹ��
void sort_unvisited_nei_kprob(vector<double>& kprob, vector<int>& index) {
    std::cout << "��ʼ sort_unvisited_nei_kprob��" << endl;
    /*auto compare = [&kprob](int i, int j) {
        return std::fabs(kprob[i] - kprob[j]) > EPSILON;
        };*/
    auto compare = [&kprob](int i, int j) -> bool {
        return std::fabs(kprob[i] - kprob[j]) > EPSILON ? true : false;
        };
    std::sort(index.begin(), index.end(), compare);

    //���������
    for (int i = 0; i < index.size(); i++) {
        std::cout << index[i] << ": " << kprob[index[i]] << endl;
    }
}


//ɾ���߼���low
double Uncertain_Core::delete_compute_low(vector<double>& kthres, int k, int root) {
    double low;
    double upper = kthres[root];
    int d = deg[root];  //ֱ���ھ����ж�
    if (d < k) {
        low = 0;
    }
    else
    {
        int count = 0;  //��¼kprob��С��root���ھ�����
        vec_b visited(d, false);
        for (int t = 0; t < d; t++) {
            int w = adj[root][t].u;
            double th = kthres[w];
            if ((th - upper) >= -EPSILON) {
                visited[t] = true;
                count++;
            }
        }

        if (count >= k) {
            low = kprob_comp_scale(root, visited, k);
        }
        else
        {
            low = 0;
        }

    }
    //cout << "root:" << root << " low:" << low << " upper:" << upper << endl;
    return low;
}


int Uncertain_Core::delete_search_candi(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiNode, vec_d& candiProb, vec_i& candiRIndex, vec_i& neiNum) {
    int count = 0;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //��ʼ����Ҫ���ʵĽڵ�
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) > EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;
    vec_b visitedCandi(n, false);
    que.push(root[0]);
    visitedCandi[root[0]] = true;
    if (root[1] != -1) {    //����ڶ���root
        que.push(root[1]);
        visitedCandi[root[1]] = true;
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        candiNode[count] = node;
        candiProb[count] = kprob_comp(node, visited, k + 1);    //��ʼ�ھ�������С��k��������=0       
        candiRIndex[node] = count;
        //cout << "candiNode:" << node << " prob:" << candiProb[count] << endl;
        count++;        

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if ((kthres[w] - low) > EPSILON) {  //kthres[w] > low
                neiNum[node]++;
                if (!visitedCandi[w] && (kthres[w] - upper) < EPSILON) {    //kthres[w] <= upper
                    que.push(w);
                    visitedCandi[w] = true;
                    //cout << "nei:" << w << endl;
                }
            }
        }
    }

    return count;
}

//���º�ѡ����thres���ҳ�ʼֵ����Ϊlow
void Uncertain_Core::delete_update_candidate_thres(vector<double>& kthres, int& k, int& count, double low, vec_i& candiNode, vec_d& candiProb, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum) {
    int index;
    double curThres = low;
    int flag = count;
    while (flag) {
        int minIndex = 0;
        int minV = candiNode[0];
        double p = candiProb[0];
        for (int i = 1; i < count; i++) {
            if (p > candiProb[i]) {
                minIndex = i;
                p = candiProb[i];
                minV = candiNode[i];
            }
        }
        curThres = max(curThres, p);
        //cout << "v:" << minV << " curThres:" << curThres << endl;
        //cout << "ԭ��v��thres��" << kthres[minV] << endl;
        kthres[minV] = curThres;
        visited[minV] = false;
        candiProb[minIndex] = 2;   //����ֵ�������2
        flag--;
        //�������ھ�
        for (int t = 0; t < deg[minV]; t++) {
            int w = adj[minV][t].u;
            index = candiRIndex[w];
            //w��candiSet�в���Ҫ����
            if (index != -1 && candiProb[index] != 2) {
                neiNum[w]--;
                if (neiNum[w] < k + 1) {
                    candiProb[index] = 0;
                }
                else
                {
                    candiProb[index] = kprob_comp(w, visited, k + 1);
                }
                //cout << "nei:" << w << " kprob:" << candiProb[index] << endl;
            }
        }
    }
}


//ɾ���ߺ� ���ݺ�ѡ�����и���
vector<vector<double> > Uncertain_Core::delete_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    std::cout << "ɾ���� & ����ԭʼͼ�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;     //��ѡ����������
    vec_i candiNode(n);
    vec_d candiProb(n);
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = delete_compute_low(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = delete_compute_low(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = delete_compute_low(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = delete_compute_low(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (k == kk - 1) {
            break;
        }
        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, candiNode, candiProb, node_visited, candiReverseIndex, node_visited_num);
        }        
    }
    
    //k = kk-1��Ѱ�Һ�ѡ�� + ����thres
    //cout << "kk-1 -" << " low:" << low << " up:" << up << endl;
    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        time = omp_get_wtime() - tm;
        return range;
    }

    //Ѱ�Һ�ѡ��
    count = delete_search_candi(thres[kk - 1], low, up, root, kk - 1, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
    ct += count;
    std::cout << "kk-count:" << count << endl;

    int index;
    int flag = count;
    double curThres = low;
    bool eraseflag = false;  //ά��kmax
    while (flag) {
        int minIndex = 0;
        int minV = candiNode[0];
        double p = candiProb[0];
        for (int i = 1; i < count; i++) {
            if (p > candiProb[i]) {
                p = candiProb[i];
                minV = candiNode[i];
                minIndex = i;
            }
        }
        curThres = max(curThres, p);
        thres[kk - 1][minV] = curThres;
        //cout << "minV:" << minV << " curThres:" << curThres << endl;
        if (fabs(curThres - 0) < EPSILON) {     //ԭ��thres��Ϊ0
            core[minV]--;
            eraseflag = true;   //kmax������Ҫ����1
        }
        node_visited[minV] = false;
        candiProb[minIndex] = 2;
        flag--;

        //�������ھ�
        for (int t = 0; t < deg[minV]; t++) {
            int w = adj[minV][t].u;
            index = candiReverseIndex[w];
            if (index != -1 && candiProb[index] != 2) {
                node_visited_num[w]--;
                if (node_visited_num[w] < kk) {
                    candiProb[index] = 0;
                }
                else
                {
                    candiProb[index] = kprob_comp(w, node_visited, kk);
                }
            }
        }
    }

    //����kmax -- ��������Ҫ��
    if (kk == kmax && eraseflag) {
        std::cout << "��ʼ��ѯkmax�Ƿ���Ҫ�ı�!!!" << endl;
        for (int t = 0; t < n; t++) {
            if (core[t] == kmax) {
                eraseflag = false;
                std::cout << "��ѯ��������" << t << endl;
                break;
            }
        }
        if (eraseflag) {
            thres.erase(thres.begin() + kmax - 1);
            kmax--;
            std::cout << "kk == kmax && kmax--" << endl;
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}



//delete�Ż�������0)��ʼ������1����̬����range��2�����ݷ������½磻3������ȷ����ѡ���thres��4���Ż������Ľ��
vector<vector<double> > Uncertain_Core::delete_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "candidate--kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));    //������range��Ϊ����ֵ
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiNode(n);
    vec_d candiProb(n);
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ -- ����
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = delete_compute_low(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = delete_compute_low(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = delete_compute_low(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = delete_compute_low(thres[k], k + 1, u);
        }
        //cout << "low:" << low << " upper:" << up << endl;
        /*range[0][k] = low;
        range[1][k] = up;*/

        if (k == kk - 1) {
            break;
        }
        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //��̬����range��Ѱ�Һ�ѡ��
        count = d_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, candiNode, candiProb, node_visited, candiReverseIndex, node_visited_num);
        }
    }

    //k = kk-1��Ѱ�Һ�ѡ�� + ����thres
    //cout << "kk-1 -" << " low:" << low << " up:" << up << endl;
    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        time = omp_get_wtime() - tm;
        return range;
    }

    //��̬����range��Ѱ�Һ�ѡ��
    count = d_candi_dynamic_update_range(thres[kk - 1], low, up, root, kk - 1, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
    ct += count;

    if (count != 0) {
        int index;
        int flag = count;
        double curThres = low;
        bool eraseflag = false;  //ά��kmax
        while (flag) {
            int minIndex = 0;
            int minV = candiNode[0];
            double p = candiProb[0];
            for (int i = 1; i < count; i++) {
                if (p > candiProb[i]) {
                    p = candiProb[i];
                    minV = candiNode[i];
                    minIndex = i;
                }
            }
            curThres = max(curThres, p);
            thres[kk - 1][minV] = curThres;
            //cout << "minV:" << minV << " curThres:" << curThres << endl;

            if (fabs(curThres - 0) < EPSILON) {     //ԭ��thres��Ϊ0
                core[minV]--;
                eraseflag = true;   //kmax������Ҫ����1
            }
            node_visited[minV] = false;
            candiProb[minIndex] = 2;
            flag--;

            //�������ھ�
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                index = candiReverseIndex[w];
                if (index != -1 && candiProb[index] != 2) {
                    node_visited_num[w]--;
                    if (node_visited_num[w] < kk) {
                        candiProb[index] = 0;
                    }
                    else
                    {
                        candiProb[index] = kprob_comp(w, node_visited, kk);
                    }
                }
            }
        }

        //����kmax -- ��������Ҫ��
        if (kk == kmax && eraseflag) {
            for (int t = 0; t < n; t++) {
                if (thres[kk - 1][t] >= EPSILON) {
                    eraseflag = false;
                    break;
                }
            }
            if (eraseflag) {
                thres.erase(thres.begin() + kmax - 1);
                kmax--;
                std::cout << "kk == kmax && kmax--" << endl;
            }
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


//Ѱ�Һ�ѡ��ʱ ��̬�����߽�
int Uncertain_Core::d_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiNode, vec_d& candiProb, vec_i& candiRIndex, vec_i& neiNum) {
    int count = 0;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //��ʼ����Ҫ���ʵĽڵ�
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) > EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;
    vec_b visitedCandi(n, false);
    que.push(root[0]);
    visitedCandi[root[0]] = true;
    if (root[1] != -1) {
        que.push(root[1]);
        visitedCandi[root[1]] = true;
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        candiNode[count] = node;
        candiProb[count] = kprob_comp(node, visited, k + 1);
        candiRIndex[node] = count;
        count++;

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if ((kthres[w] - low) > EPSILON) {
                neiNum[node]++;
                if (!visitedCandi[w] && (kthres[w] - kthres[node]) < EPSILON) {     //��������upper
                    que.push(w);
                    visitedCandi[w] = true;
                }
            }
        }
    }

    return count;
}


//ɾ���� �Ƚ��������������low
void Uncertain_Core::delete_compare_range(vector<vector<double> >& thres, double scale) {
    double tm_1;
    double tm_2;
    int count1, count2;    
    selectRandomEdges(scale);

    int um = unselected.size();
    for (int i = 0; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        std::cout << "��i���ߣ�" << i << " u:" << u << " v:" << v << endl;

        //ά��ͼ�ṹ
        bool flag = false;
        int d = deg[u];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[u][j].u;
                if (w == v) {
                    flag = true;
                }
            }
            else
            {
                adj[u][j - 1].u = adj[u][j].u;
                adj[u][j - 1].p = adj[u][j].p;
            }
        }

        flag = false;
        d = deg[v];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[v][j].u;
                if (w == u) {
                    flag = true;
                }
            }
            else
            {
                adj[v][j - 1].u = adj[v][j].u;
                adj[v][j - 1].p = adj[v][j].p;
            }
        }
        deg[u]--;
        deg[v]--;


        int mincore = min(core[u], core[v]);
        //��ʼ�����½�
        vector<vector<double> > range1 = delete_threshold_compute_candidate(thres, u, v, tm_1, count1);
        //���ݷ���ȡ�½�
        vector<vector<double> > range2 = delete_threshold_compute_shrink_low(thres, u, v, tm_1, count2);

        delete[] core;
		vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);

        //�Ƚ�ɾ���ߺ�ı仯��ֵ
        compareArrays(thres, thres_2, range2[0], range2[1], mincore - 1);   //���ܱȽ�core��С�ĵ�

        for (int i = 0; i < range1[0].size(); i++) {
            if ((range2[0][i] - range1[0][i]) > EPSILON) {
                std::cout << "low���� ------- k��" << i;
                std::cout << "  initial��" << range1[0][i] << " opt:" << range2[0][i] << endl;
            }
        }

        //��newArray���Ƹ�thres
        if (thres.size() == thres_2.size()) {
            thres = thres_2;
        }
        else
        {
            thres.resize(thres_2.size());
            for (size_t i = 0; i < kmax; i++) {
                thres[i] = thres_2[i];
            }
        }
    }
    unselected.clear();
}


//�Ż����� -- ���ݷ���ȡ�½�
vector<vector<double> > Uncertain_Core::delete_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "candidate--kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));    //������range��Ϊ����ֵ
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiNode(n);
    vec_d candiProb(n);
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root

        //���ݷ������½�
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = d_shrink_lower_bound(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = d_shrink_lower_bound(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = d_shrink_lower_bound(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = d_shrink_lower_bound(thres[k], k + 1, u);
        }
        //cout << "low:" << low << " upper:" << up << endl;
        range[0][k] = low;
        range[1][k] = up;

        if (k == kk - 1) {
            break;
        }
        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
        ct += count;

        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, candiNode, candiProb, node_visited, candiReverseIndex, node_visited_num);
        }
    }

    //k = kk-1��Ѱ�Һ�ѡ�� + ����thres
    //cout << "kk-1 -" << " low:" << low << " up:" << up << endl;
    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        time = omp_get_wtime() - tm;
        return range;
    }

    //Ѱ�Һ�ѡ��
    count = delete_search_candi(thres[kk - 1], low, up, root, kk - 1, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
    ct += count;

    if (count != 0) {
        int index;
        int flag = count;
        double curThres = low;
        bool eraseflag = false;  //ά��kmax
        while (flag) {
            int minIndex = 0;
            int minV = candiNode[0];
            double p = candiProb[0];
            for (int i = 1; i < count; i++) {
                if (p > candiProb[i]) {
                    p = candiProb[i];
                    minV = candiNode[i];
                    minIndex = i;
                }
            }
            curThres = max(curThres, p);
            thres[kk - 1][minV] = curThres;
            //cout << "minV:" << minV << " curThres:" << curThres << endl;

            if (fabs(curThres - 0) < EPSILON) {     //ԭ��thres��Ϊ0
                core[minV]--;
                eraseflag = true;   //kmax������Ҫ����1
            }
            node_visited[minV] = false;
            candiProb[minIndex] = 2;
            flag--;

            //�������ھ�
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                index = candiReverseIndex[w];
                if (index != -1 && candiProb[index] != 2) {
                    node_visited_num[w]--;
                    if (node_visited_num[w] < kk) {
                        candiProb[index] = 0;
                    }
                    else
                    {
                        candiProb[index] = kprob_comp(w, node_visited, kk);
                    }
                }
            }
        }

        //����kmax -- ��������Ҫ��
        if (kk == kmax && eraseflag) {
            for (int t = 0; t < n; t++) {
                if (thres[kk - 1][t] >= EPSILON) {
                    eraseflag = false;
                    break;
                }
            }
            if (eraseflag) {
                thres.erase(thres.begin() + kmax - 1);
                kmax--;
                std::cout << "kk == kmax && kmax--" << endl;
            }
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


//��һ�������½� -- ��Ϊopt1��opt2 ���ߵ���ʵ�� �Ż����ܲ�ͬ
double Uncertain_Core::d_shrink_lower_bound(vector<double>& kthres, int k, int root) {
    double low;
    double upper = kthres[root];
    int d = deg[root];  //ֱ���ھ����ж�
    if (d < k) {
        low = 0;
    }
    else
    {
        int coreD = 0;  //����core��С��k���ھ�����
        int count = 0;  //��¼kprob��С��root���ھ�����
        vec_b visited(d, false);
        vector<pair<double, int> > unvisitedNei; //��¼kprobС��root���ھ�
        for (int t = 0; t < d; t++) {
            int w = adj[root][t].u;
            double th = kthres[w];
            if ((th - upper) >= -EPSILON) {     //th >= EPSILON
                visited[t] = true;
                coreD++;
                count++;
            }
            else
            {
                if (core[w] >= k) {
                    unvisitedNei.push_back(make_pair(th, t));
                    coreD++;
                }
            }
        }

        //d >= coreD >= count
        if (coreD < k) {
            low = 0;
        }
        else if (count >= k) {  //�����ھ������½�
            low = kprob_comp_scale(root, visited, k);

            if (coreD > count) {    //ȷ��unvisitedNei�����ھ� -- opt1
                sort(unvisitedNei.begin(), unvisitedNei.end(), [](const pair<double, int>& p1, const pair<double, int>& p2) {
                    return p1.first > p2.first;
                    });

                double th = unvisitedNei[0].first;
                if ((th - low) > EPSILON) {
                    bool flag = true;
                    int countNum = count;
                    while (flag) {
                        int index = unvisitedNei[countNum - count].second;
                        int w = adj[root][index].u;
                        visited[index] = true;
                        double kprob = kprob_comp_scale(root, visited, k);
                        countNum++;
                        if ((kprob - kthres[w]) >= -EPSILON) {    //ȷ��lowȡֵ
                            low = kthres[w];
                            flag = false;
                        }
                        else
                        {
                            low = kprob;    //�ܷ�һֱ��ֵ������
                            if ((countNum == coreD) || ((kprob - unvisitedNei[countNum - count].first) >= -EPSILON)) {      //��ζ������һ���� kprob��Ȼ����
                                flag = false;
                            }
                        }
                    }
                }
            }            
            
        }
        else
        {
            //����ѡ��ʣ�µ��ھ� -- opt2
            sort(unvisitedNei.begin(), unvisitedNei.end(), [](const pair<double, int>& p1, const pair<double, int>& p2) {
                return p1.first > p2.first;
                });

            bool flag = true;
            int countNum = count;
            while (flag) {
                int index = unvisitedNei[countNum - count].second;
                int w = adj[root][index].u;
                visited[index] = true;
                double kprob = kprob_comp_scale(root, visited, k);
                countNum++;
                if ((kprob - kthres[w]) >= -EPSILON) {    //ȷ��lowȡֵ
                    low = kthres[w];
                    flag = false;
                }
                else
                {
                    low = 0;
                    if (countNum == coreD) {
                        flag = false;
                    }
                }
            }
        }
    }
    //cout << "root:" << root << " low:" << low << " upper:" << upper << endl;
    return low;
}


//�������º�ѡ��
vector<vector<double> > Uncertain_Core::delete_threshold_compute_batchUp(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    vec_i candiNode(n);
    vec_d candiProb(n);
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = delete_compute_low(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = delete_compute_low(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = delete_compute_low(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = delete_compute_low(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (k == kk - 1) {
            break;
        }
        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            //delete_update_candidate_thres(thres[k], k, count, low, candiNode, candiProb, node_visited, candiReverseIndex, node_visited_num);
            d_batch_update_candidate(thres[k], k, count, low, candiNode, candiProb, node_visited, candiReverseIndex, node_visited_num);
        }
    }

    //k = kk-1��Ѱ�Һ�ѡ�� + ����thres
    //cout << "kk-1 -" << " low:" << low << " up:" << up << endl;
    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        time = omp_get_wtime() - tm;
        return range;
    }

    //Ѱ�Һ�ѡ��
    count = delete_search_candi(thres[kk - 1], low, up, root, kk - 1, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
    std::cout << "kk-count:" << count << endl;

    std::queue<int> que;
    int index;
    int flag = count;
    double curThres = low;
    bool eraseflag = false;  //ά��kmax
    while (flag) {
        int minIndex, minV;
        double p = 2;
        for (int i = 0; i < count; i++) {     
            if ((candiProb[i] - curThres) < EPSILON) {
                que.push(i);
                //std::cout << "que+:" << i << " candiProb:" << candiProb[i] << endl;
            }
            if (que.empty()) { //û��С��curThres�ĵ� -- �ҵ�һ����С��
                if (p > candiProb[i]) {
                    minIndex = i;
                    p = candiProb[i];
                    minV = candiNode[i];
                }
            }            
        }

        if (!que.empty()) {
            std::unordered_set<int> neiSet;
            while (!que.empty()) {
                index = que.front();
                que.pop();

                int node = candiNode[index];
                thres[kk - 1][node] = curThres;
                //std::cout << "�������£�" << node << " thres:" << curThres << endl;
                if (fabs(curThres - 0) < EPSILON) {     //ԭ��thres��Ϊ0
                    core[node]--;
                    eraseflag = true;   //kmax������Ҫ����1
                }
                node_visited[node] = false;
                candiProb[index] = 2;
                flag--;

                //��Ҫ���µ��ھ�
                for (int t = 0; t < deg[node]; t++) {
                    int w = adj[node][t].u;
                    index = candiReverseIndex[w];
                    //w��candiSet�в���Ҫ����
                    if (index != -1 && candiProb[index] != 2) {
                        node_visited_num[w]--;
                        neiSet.insert(w);
                    }
                }
            }

            for (auto it = neiSet.begin(); it != neiSet.end(); it++) {
                int w = *it;
                if (node_visited[w]) {
                    index = candiReverseIndex[w];
                    if (node_visited_num[w] < kk) {
                        candiProb[index] = 0;
                    }
                    else
                    {
                        candiProb[index] = kprob_comp(w, node_visited, kk);
                    }
                    //std::cout << "���������ھӣ�" << w << " candiProb:" << candiProb[index] << endl;
                }                
            }
        }
        else
        {
            curThres = max(curThres, p);
            thres[kk - 1][minV] = curThres;
            //cout << "minV:" << minV << " curThres:" << curThres << endl;
            if (fabs(curThres - 0) < EPSILON) {     //ԭ��thres��Ϊ0
                core[minV]--;
                eraseflag = true;   //kmax������Ҫ����1
            }
            node_visited[minV] = false;
            candiProb[minIndex] = 2;
            flag--;

            //�������ھ�
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                index = candiReverseIndex[w];
                if (index != -1 && candiProb[index] != 2) {
                    node_visited_num[w]--;
                    if (node_visited_num[w] < kk) {
                        candiProb[index] = 0;
                    }
                    else
                    {
                        candiProb[index] = kprob_comp(w, node_visited, kk);
                    }
                }
            }
        }        
    }

    //����kmax -- ��������Ҫ��
    if (kk == kmax && eraseflag) {
        std::cout << "��ʼ��ѯkmax�Ƿ���Ҫ�ı�!!!" << endl;
        for (int t = 0; t < n; t++) {
            if (core[t] == kmax) {
                eraseflag = false;
                std::cout << "��ѯ��������" << t << endl;
                break;
            }
        }
        if (eraseflag) {
            thres.erase(thres.begin() + kmax - 1);
            kmax--;
            std::cout << "kk == kmax && kmax--" << endl;
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


//�������º�ѡ��
void Uncertain_Core::d_batch_update_candidate(vector<double>& kthres, int& k, int& count, double low, vec_i& candiNode, vec_d& candiProb, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum) {
    std::queue<int> que;    //����que�������� p < curThres �ĵ�
    int index;
    double curThres = low;
    int flag = count;
    while (flag) {
        int minIndex, minV;
        double p = 2;
        for (int i = 0; i < count; i++) {  
            if ((candiProb[i] - curThres) < EPSILON) {
                que.push(i);
            }
            if (que.empty()) { //û��С��curThres�ĵ� -- �ҵ�һ����С��
                if (p > candiProb[i]) {
                    minIndex = i;
                    p = candiProb[i];
                    minV = candiNode[i];
                }
            }                                    
        }

        if (!que.empty()) { //thres��ΪcurThres
            std::unordered_set<int> neiSet; //�����¶���
            while (!que.empty()) {
                index = que.front();
                que.pop();

                int node = candiNode[index];
                kthres[node] = curThres;
                visited[node] = false;
                candiProb[index] = 2;
                flag--;

                //��Ҫ���µ��ھ�
                for (int t = 0; t < deg[node]; t++) {
                    int w = adj[node][t].u;
                    index = candiRIndex[w];
                    //w��candiSet�в���Ҫ����
                    if (index != -1 && candiProb[index] != 2) {
                        neiNum[w]--;
                        neiSet.insert(w);
                    }
                }
            }

            for (auto it = neiSet.begin(); it != neiSet.end(); it++) {
                int w = *it;
                if (visited[w]) {
                    index = candiRIndex[w];
                    //std::cout << "���������ھӣ�" << w << endl;
                    if (neiNum[w] < k + 1) {
                        candiProb[index] = 0;
                    }
                    else
                    {
                        candiProb[index] = kprob_comp(w, visited, k + 1);
                    }
                }                
            }
        }
        else
        {
            curThres = max(curThres, p);

            kthres[minV] = curThres;
            visited[minV] = false;
            candiProb[minIndex] = 2;   //����ֵ�������2
            flag--;
            //�������ھ�
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                index = candiRIndex[w];
                //w��candiSet�в���Ҫ����
                if (index != -1 && candiProb[index] != 2) {
                    neiNum[w]--;
                    if (neiNum[w] < k + 1) {
                        candiProb[index] = 0;
                    }
                    else
                    {
                        candiProb[index] = kprob_comp(w, visited, k + 1);
                    }
                    //cout << "nei:" << w << " kprob:" << candiProb[index] << endl;
                }
            }
        }
    }
}


vector<vector<double> > Uncertain_Core::delete_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    std::cout << "ɾ���� & ����ԭʼͼopt�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiNode(n);
    vec_d candiProb(n);
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = d_shrink_lower_bound(thres[k], k + 1, u);   //����ʵ��Ч�� ����ѡ��һ��������ʽ

            root[1] = v;
            double kprob2 = d_shrink_lower_bound(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = d_shrink_lower_bound(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = d_shrink_lower_bound(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (k == kk - 1) {
            break;
        }
        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //��̬Ѱ�Һ�ѡ��
        count = d_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, candiNode, candiProb, node_visited, candiReverseIndex, node_visited_num);    //�Ƿ��������£�����
        }
    }

    //k = kk-1��Ѱ�Һ�ѡ�� + ����thres
    //cout << "kk-1 -" << " low:" << low << " up:" << up << endl;
    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        time = omp_get_wtime() - tm;
        return range;
    }

    //Ѱ�Һ�ѡ��
    count = d_candi_dynamic_update_range(thres[kk - 1], low, up, root, kk - 1, node_visited, candiNode, candiProb, candiReverseIndex, node_visited_num);
    ct += count;
    std::cout << "kk-count:" << count << endl;

    //�Ƿ��������£�����
    int index;
    int flag = count;
    double curThres = low;
    bool eraseflag = false;  //ά��kmax
    while (flag) {
        int minIndex = 0;
        int minV = candiNode[0];
        double p = candiProb[0];
        for (int i = 1; i < count; i++) {
            if (p > candiProb[i]) {
                p = candiProb[i];
                minV = candiNode[i];
                minIndex = i;
            }
        }
        curThres = max(curThres, p);
        thres[kk - 1][minV] = curThres;
        //cout << "minV:" << minV << " curThres:" << curThres << endl;
        if (fabs(curThres - 0) < EPSILON) {     //ԭ��thres��Ϊ0
            core[minV]--;
            eraseflag = true;   //kmax������Ҫ����1
        }
        node_visited[minV] = false;
        candiProb[minIndex] = 2;
        flag--;

        //�������ھ�
        for (int t = 0; t < deg[minV]; t++) {
            int w = adj[minV][t].u;
            index = candiReverseIndex[w];
            if (index != -1 && candiProb[index] != 2) {
                node_visited_num[w]--;
                if (node_visited_num[w] < kk) {
                    candiProb[index] = 0;
                }
                else
                {
                    candiProb[index] = kprob_comp(w, node_visited, kk);
                }
            }
        }
    }

    //����kmax -- ��������Ҫ��
    if (kk == kmax && eraseflag) {
        std::cout << "��ʼ��ѯkmax�Ƿ���Ҫ�ı�!!!" << endl;
        for (int t = 0; t < n; t++) {
            if (core[t] == kmax) {
                eraseflag = false;
                std::cout << "��ѯ��������" << t << endl;
                break;
            }
        }
        if (eraseflag) {
            thres.erase(thres.begin() + kmax - 1);
            kmax--;
            std::cout << "kk == kmax && kmax--" << endl;
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}

