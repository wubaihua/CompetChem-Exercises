#include <iostream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
//typedef Matrix VectorXi;
using namespace std;
//using namespace Eigen;
void showmatrix(Matrix);

int main()
{
    ifstream geo("geom.dat");
    int natom;
    geo >> natom;

    double *zval=new double[natom];
    double *x=new double[natom];
    double *y=new double[natom];
    double *z=new double[natom];

    for(int i=0;i<natom;i++)
        geo >> zval[i]>>x[i]>>y[i]>>z[i];

    geo.close();

    cout<<"number of atoms="<<natom<<endl;
    cout<<"Z,x,y,z:\n";
    for(int i=0;i<natom;i++)
        {printf("%20.12f %20.12f %20.12f %20.12f\n",zval[i],x[i],y[i],z[i]);}
    
    int nele;
    nele=0;
    for(int i=0;i<natom;i++)
    {
        nele+=zval[i];
    }
    cout<<"nele="<<nele<<endl;


    double* bond = new double[natom,natom];
    for(int i=0;i<natom;i++)
    {
        for(int j=0;j<natom;j++)
        {
            bond[i,j]=sqrt(pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2)+pow((z[i]-z[j]),2));
            //cout<<i<<"-"<<j<<"="<<bond[i,j];
        }     
    }
        
    /*
    ifstream hessian("hess.dat");
    hessian >> natom;
    Matrix hess(3*natom,3*natom);
    for(int i=0;i<natom*3;i++)
    {
        for(int j=0;j<natom;j++)
        {
            hessian >> hess(i,3*j)>>hess(i,3*j+1)>>hess(i,3*j+2);
            cout<<"hess"<<hess(i,3*j)<<hess(i,3*j+1)<<hess(i,3*j+2)<<endl;
        }
    }
    hessian.close();

    double *mass=new double[3*natom];
    Matrix whess(3*natom,3*natom);
    for(int i=0;i<natom;i++)
    {
        mass[3*i]=zval[i];
        mass[3*i+1]=zval[i];
        mass[3*i+2]=zval[i];
        cout<<"mass"<<mass[3*i]<<mass[3*i+1]<<mass[3*i+2]<<endl;
    }

    for(int i=0;i<3*natom;i++)
    {
        for(int j=0;j<3*natom;j++)
        {
            whess(i,j)=hess(i,j)/sqrt(mass[i]*mass[j]);
            cout<<i<<"-"<<j<<"="<<whess(i,j)<<endl;
        }
    }


    
    Eigen::SelfAdjointEigenSolver<Matrix> solver(whess);
    Matrix evecs=solver.eigenvectors();
    Matrix evals=solver.eigenvalues();

    for(int i=0;i<3*natom;i++)
    {
        cout<<evals(i)*sqrt(4.35873E35/(pow(5.291772,2)*9.10938))/(2*3.14159265*3.0E8*100)<<endl;
    }
    */



    










    ifstream np("enuc.dat");
    double npe;
        np >> npe;

    np.close();




    ifstream s1("s.dat");
    int nbasis;
    int q;double a;
    bool iq;
    while (!s1.eof())
    {
        s1 >> nbasis >> q >> a;
        //cout << nbasis << endl;
    }
    s1.close();
    cout << nbasis << endl;

    s1.close();
    /*

    char c;
    int nline=0;
    int nbasis;
    ifstream s1("s.dat");
    string c2;
  `
    while(getline(s1,c2))
    {
        /*if(s1.peek()==EOF)
        {
            break;
        }
        //if(c=='\n');
        //cout<<s1.eof(); 
        s1 >> c2;
        nbasis=nbasis+1;
    }
    
    s1.close();
    */
    ifstream s2("s.dat");
    int k;int l;double p;
    for(int i=0;i<nbasis;i++)
    {
        s2 >>k>>l>>p;
        //cout<<k<<l<<p<<endl;
    }
    //nbasis=k;
    cout<<"nbasis="<<nbasis<<endl;
    s2.close();
    
    Matrix s(nbasis,nbasis);
    ifstream s3("s.dat");
    for(int i=0;i<nbasis;i++)
    {
        for(int j=0;j<=i;j++)
        {
            s3 >> k >> l >> s(i,j);
            //cout<<i<<"="<<j<<"="<<s(i,j)<<endl;
        }
    }
    s3.close();

    

    Matrix t(nbasis,nbasis);
    ifstream t1("t.dat");
    for(int i=0;i<nbasis;i++)
    {
        for(int j=0;j<=i;j++)
        {
            t1 >> k >> l >> t(i,j);
            //cout<<s(i,j)<<endl;
        }
    }
    t1.close();

    Matrix v(nbasis,nbasis);
    ifstream v1("v.dat");
    for(int i=0;i<nbasis;i++)
    {
        for(int j=0;j<=i;j++)
        {
            v1>> k >> l >> v(i,j);
            //cout<<s(i,j)<<endl;
        }
    }
    v1.close();
    
    for(int i=0;i<nbasis;i++)
    {
        for(int j=i;j<nbasis;j++)
        {
            s(i,j)=s(j,i);
            v(i,j)=v(j,i);
            t(i,j)=t(j,i);
        }
    }
    //cout<<t(5,6)<<endl;

    Matrix Hc(nbasis,nbasis);
    Hc=v+t;

    
    ifstream di("eri.dat");
    double *eri=new double[nbasis,nbasis,nbasis,nbasis]();
    double din;
    int i;int j;
    while (!di.eof())
    {
        di >> i >> j >> k >> l >> din;
        i=i-1;j=j-1;k=k-1;l=l-1;
        eri[i,j,k,l]=din;
        //cout<<i<<"-"<<j<<"-"<<k<<"-"<<l<<"="<<eri[i,j,k,l]<<endl;
        eri[i,j,l,k]=eri[i,j,k,l];
        eri[j,i,l,k]=eri[i,j,k,l];
        eri[k,l,i,j]=eri[i,j,k,l];
        eri[k,l,j,i]=eri[i,j,k,l];
        eri[l,k,i,j]=eri[i,j,k,l];
        eri[l,k,j,i]=eri[i,j,k,l];

        
    }
    
    //cout<<eri[5,5,5,5];
    di.close();


    Eigen::SelfAdjointEigenSolver<Matrix> solver(s);
    Matrix ls=solver.eigenvectors();
    Matrix lambdas=solver.eigenvalues();
    
    //Matrix ls=s.Eigenvectors();
    //Matrix lambdas=s.eigenvalues();
    

    Matrix lambda_m=Matrix::Zero(nbasis,nbasis);
    

    for (int i = 0; i <nbasis; i++)
    {
        lambda_m(i,i)=pow(lambdas(i),-0.5);
    }
    
    Matrix s_m(nbasis,nbasis);//S^{-1/2}
    s_m=ls*lambda_m*ls.adjoint();
    /**
    for (int i = 0; i < nbasis; i++)
    {
        cout<<s_m.row(i)<<endl;
    }
    */

    Matrix F0=s_m.adjoint()*Hc*s_m;
    /*
    for (int i = 0; i < nbasis; i++)
    {
        cout<<F0.row(i)<<endl;
    }
    */
    Eigen::SelfAdjointEigenSolver<Matrix> g(F0);
    Matrix C0=g.eigenvectors();
    Matrix eps=g.eigenvalues();
    
    Matrix C=s_m*C0;

    //showmatrix(C);

    Matrix D=Matrix::Zero(nbasis,nbasis);
    double sum1;
    for (int mu = 0; mu < nbasis; mu++)
    {
        for (int nu = 0; nu < nbasis; nu++)
        {
            sum1=0;
            for (int m = 0; m < nele/2; m++)
            {
                sum1+=C(mu,m)*C(nu,m);
            }
            D(mu,nu)=sum1;
        }
        
    }
    //showmatrix(F0);
    //showmatrix(D);
    //cout<<"nele="<<nele<<endl;
    double E_ele=0;
    for (int mu = 0; mu < nbasis; mu++)
    {
        for (int nu = 0; nu < nbasis; nu++)
        {
            E_ele=E_ele+D(mu,nu)*(Hc(mu,nu)+F0(mu,nu));
        }
        
    }
    
    cout<<"E_ele_0="<<E_ele<<endl;

    double E_t=E_ele+npe;

    cout<<"E_total_0="<<E_t<<endl;

    double deltaE=100;
    double rmsd=100;
    int istep=1;
    double E_ele0;
    Matrix D0;

    while (deltaE>1.0E-6 || rmsd>1.0E-6)
    {
        E_ele0=E_ele;
        D0=D;
        for (int mu = 0; mu < nbasis; mu++)
        {   
            for (int nu = 0; nu < nbasis; nu++)
            {
                double sum2=0;
                for (int lambda = 0; lambda < nbasis; lambda++)
                {
                    for (int sigma = 0; sigma < nbasis; sigma++)
                    {       
                        sum2=sum2+D0(lambda,sigma)*(2*eri[mu,nu,lambda,sigma]-eri[mu,lambda,nu,sigma]);

                    }
                
                }
                F0(mu,nu)=Hc(mu,nu)+sum2;
            }
        
        }

        Matrix F=s_m.adjoint()*F0*s_m;

        Eigen::SelfAdjointEigenSolver<Matrix> g1(F);
        Matrix Cp=g1.eigenvectors();
        eps=g1.eigenvalues();

        C=s_m*Cp;

        for (int mu = 0; mu < nbasis; mu++)
        {
            for (int nu = 0; nu < nbasis; nu++)
            {
                double sum3=0;
                for (int m = 0; m < nele/2; m++)
                {
                    sum3+=C(mu,m)*C(nu,m);
                }
                D(mu,nu)=sum3;
            }
            
        }
        E_ele=0;
       
        for (int mu = 0; mu < nbasis; mu++)
        {
            for (int nu = 0;nu < nbasis; nu++)
            {
                E_ele+=D(nu,mu)*(Hc(mu,nu)+F(mu,nu));
            }
            
        }
        
    //    for (int i = 0; i < nele/2; i++)
    //    {
    //        E_ele+=2*eps(i);
    //    }
       
        
        deltaE=abs(E_ele-E_ele0);
        rmsd=0;
        for (int mu = 0; mu < nbasis; mu++)
        {
            for (int nu = 0; nu < nbasis; nu++)
            {
               rmsd+=pow(D(mu,nu)-D0(mu,nu),2);
            }
            
        }
        rmsd=pow(rmsd,0.5);

        cout<<"step="<<istep<<endl;
        cout<<"E_ele="<<E_ele<<endl;
        cout<<"E_tot="<<E_ele+npe<<endl;
        cout<<"deltaE="<<deltaE<<endl;
        cout<<"RMSD="<<rmsd<<endl;

        istep+=1;


    }
    

















    return 0;

}

void showmatrix(Matrix A)
{
    int n=A.rows();
    for (int i = 0; i <n; i++)
    {
        cout<<A.row(i)<<endl;
    }
    
}
