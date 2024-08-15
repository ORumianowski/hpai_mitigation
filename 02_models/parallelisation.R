

library(Rcpp)

cppFunction('
            
    #include <iostream>
    #include <vector>
    #include <thread>        
            
    void sumPart(const std::vector<int>& data, int start, int end, int& result) {
    result = 0;
    for (int i = start; i < end; ++i) {
        result += data[i];
    }
}

int main() {
    std::vector<int> data(1000000, 1);  
int result1 = 0, result2 = 0;

// Créer deux threads pour diviser la tâche
std::thread t1(sumPart, std::ref(data), 0, data.size() / 2, std::ref(result1));
std::thread t2(sumPart, std::ref(data), data.size() / 2, data.size(), std::ref(result2));

// Attendre que les threads terminent
t1.join();
t2.join();

// Calculer le résultat final
int finalResult = result1 + result2;
std::cout << "La somme est: " << finalResult << std::endl;

return 0;
}')

