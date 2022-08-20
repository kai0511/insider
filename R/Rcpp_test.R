cppFunction("Eigen::VectorXd mult(Eigen::VectorXd v) {
    Eigen::VectorXd m = (v.array() > 0).select;
    return m;
}", depends = "RcppEigen")

cppFunction("#include <iostream>
Eigen::VectorXi v = Eigen::VectorXi::Random(4);
cout << 'Here is the vector v:\n';
for(auto x : v) cout << x << ' ';
cout << '\n';", depends = "RcppEigen")



