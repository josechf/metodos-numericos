
#include <bits/stdc++.h>

#define stop system("pause");
#define cls system("cls");

using  namespace std;
void print_matriz(float[][6],int,float[]);
void lineales();
void cuadraticas();

// CLASES -------------
class Determinant{
    public:
    float matriz[6][6],independent[6];
    int order;
    public:
    void get_data();
    float determinant(float [][6],int);
    float cofactors(float [][6],int,int,int);
    bool divergir(float[][6],int);
};

class Gauss : Determinant{
    public:
    void eliminacion(float[][6],float[],int);
    void jordan(float[][6],float[],int);
    void seidel(float [][6],float[],int);
};

class Raices{
    public:
    double rank_left,rank_right,approach[1000],tolerance; int iteration,exponencial; bool fin;
    float coefficient[10],constant;
    float coefficient_der[10],constant_der;
    public:
    Raices();
    void cal_newton(double,int);
    double newton(double);
    double derivada(double);
    double biseccion(double,double);
    double falsa(double,double);
    double calcular(char);
    double f(double);
    void get_rank_left();
    void get_rank_right();
    void get_tolerance();
    void get_polinomio();
    void get_derivada();
    int trunca(float,int);
};

    //methods Raices
Raices::Raices(){
    int i=0;
    rank_left=rank_right=tolerance=exponencial=constant=constant_der=iteration=0;fin=false;
    for(i=0;i<1000;i++){approach[i]=0;if(i<10){coefficient[i]=0;coefficient_der[i]=0;}}
}

void Raices::cal_newton(double raiz,int desplegarmenu){
    char opciones;
    iteration=0;
    if(desplegarmenu==0){
        cout <<" |A| METODO PARA ACERCARSE A LA RAIZ\n |N| NEWTON DIRECTO : "<<endl;cin>>opciones;
    }else{
        opciones='s'; //la s para meterle algo ggg
    }
    if(opciones=='a'){
        get_rank_left();
        get_rank_right();
        get_tolerance();
        get_polinomio();
        cin.clear();fflush(stdin);
        cout << " METODO\n |b|BISECCION\n |f|FALSA\n "<<endl; cin >>opciones;
        raiz = calcular(opciones);
    }
    if(opciones=='n'){
        get_rank_right();
        raiz = rank_right;
        get_polinomio();
    }
    get_derivada();
    get_tolerance();
                        // HERE BEGINS NEWTON, BEFORE HE WAS GETTING THINGS !!!!
    cout << " i\t\tA\t\tf(A)\t\tf'(A)\t\tAi+1\t\t error " << endl;
    do{
        printf(" %d\t%.11lf\t%.11lf\t%.12lf\t%.11lf\t%.11lf\n",iteration,raiz,f(raiz),derivada(raiz),newton(raiz),newton(raiz)-raiz/newton(raiz));
        if(trunca(f(raiz),7)==0){fin = true;cout << " LA RAIZ ES : "<<raiz<<endl;}
        raiz = newton(raiz);
        iteration++;
    }while(fin==false);
}

double Raices::newton(double raiz){
    return raiz - f(raiz)/derivada(raiz);
}

double Raices::derivada(double unknow){
    double suma=0;int aux=exponencial;
    for(int i=0;i<aux-1;i++){
        suma+=coefficient_der[i]*pow(unknow,exponencial-1);
        exponencial--;
    }
    exponencial=aux;
    return suma + constant_der;
    //return exp(-unknow)+1;
}

void Raices::get_derivada(){
    int x=0;
    int aux=exponencial-1;
    cout <<"\n DERIVADA"<<endl;
    for(int i=0;i<exponencial-1;i++){
        do{
            cin.clear();
            fflush(stdin);
            cout << " COEFICIENTE"<<aux<<"  : "<<endl;
            x=scanf("%f",&coefficient_der[i]);
        }while(x==0);
        aux--;
    }
    do{
        cin.clear();
        fflush(stdin);
        cout << " CONSTANTE : "<<endl;
        x=scanf("%f",&constant_der);
    }while(x==0);
}

void Raices::get_rank_right(){
    int x=0;
    do{
        cin.clear();
        fflush(stdin);
        cout << " INTERVALO DERECHO : "<<endl;
        x=scanf("%lf",&rank_right);
    }while(x==0);
}
void Raices::get_rank_left(){
int x=0;
    do{
        cin.clear();
        fflush(stdin);
        cout << " INTERVALO IZQUIERDO :  "<<endl;
        x=scanf("%lf",&rank_left);
    }while(x==0);
}
void Raices::get_polinomio(){
    int x=0;
    do{
        cin.clear();
        fflush(stdin);
        cout << " MAXIMO EXPONENTE : "<<endl;
        x=scanf("%d",&exponencial);
    }while(x==0);
    int aux=exponencial;
    for(int i=0;i<exponencial;i++){
        do{
            cin.clear();
            fflush(stdin);
            cout << " COEFICIENTE"<<aux<<"  : "<<endl;
            x=scanf("%f",&coefficient[i]);
        }while(x==0);
        aux--;
    }
    do{
        cin.clear();
        fflush(stdin);
        cout << " CONSTANTE : "<<endl;
        x=scanf("%f",&constant);
    }while(x==0);
}
void Raices::get_tolerance(){
    int x=0;
    do{
        cin.clear();
        fflush(stdin);
        cout << " TOLERANCIA : "<<endl;
        x=scanf("%lf",&tolerance);
    }while(x==0);
}

double Raices::calcular(char method){
    double iz=rank_left,der=rank_right,raiz=0;
    char opciones; int extra=0;
    extra = (method=='b') ? (log(rank_right-rank_left)-log(tolerance))/log(2) : 0;
    if(f(rank_left)*f(rank_right)>0){
        cout << " APROXIMACION INCORRECTA "<<endl;
    }else{
        if(method=='b'){ cout << " EL NUMERO DE ITERACIONES DEBE SER : "<<(log(rank_right-rank_left)-log(tolerance))/log(2)<<endl;}
        cout << "\n i\ta\t\tb\t\txm\t\tf(a)\t\tf(xm)\t\tf(a)*f(xm)\t(b-a)/b"<<endl;
        do{
            if(method=='b'){
                approach[iteration] = biseccion(rank_left,rank_right);
            }else{
                if(method=='f'){
                    approach[iteration] = falsa(rank_left,rank_right);
                }
            }
            printf(" %d\t%.11lf\t%.11lf\t%.12lf\t%.11lf\t%.11lf\t%.11lf\t%.11lf\n",iteration+1,rank_left,rank_right,approach[iteration],f(rank_left),f(approach[iteration]),f(rank_left)*f(approach[iteration]),(rank_right-rank_left)/rank_right);
            if((rank_right-rank_left)/rank_right<=tolerance){
                if(trunca(f(approach[iteration]),4)!=0){
                    cout<<" |m| MENOR TOLERANCIA : "<<endl;
                    cin >>opciones;
                    if(opciones=='m'){
                        get_tolerance();
                        if(method=='b'){ cout << " EL NUMERO DE ITERACIONES DEBE SER : "<<(log(rank_right-rank_left)-log(tolerance))/log(2)<<endl;}
                    }else{
                        fin = true;
                    }
                }
            }
            (f(rank_left)*f(approach[iteration])<0 ? rank_right : rank_left) = approach[iteration];
            /*if(method=='b'&&iteration>extra){
                if(trunca(f(approach[iteration]),4)!=0){
                    cout<<" |m| MAS ITERACIONES : "<<endl;
                    cin >>opciones;
                    if(opciones=='m'){
                        int anadir=0;
                        cout <<" ITERACIONES EXTRA : "<<endl;
                        cin >>anadir;
                        extra +=anadir;
                    }else{
                        fin = true;
                    }
                }
            }*/
            if(trunca(f(approach[iteration]),4)==0){                           //////////////POSIBLEMENTE QUITE ESTO POR NO SABER CUANTOS DECIMALES USAR
                fin = true; cout << "\n FIN DE ITERACIONES POR IGUALAR A CERO "<<endl;raiz = approach[iteration];
            }else{
               if(trunca(f(approach[iteration-1]),9)>trunca(f(approach[iteration]),9)&&trunca(f(approach[iteration]),4)>=0){ //el que mas se acerca a cero
                    raiz = approach[iteration];
               }
            }
            iteration++;
        }while(fin==false);
        if(raiz==0){
            cout <<"\nAPROXIMACION : "<<approach[iteration-1]<<endl;
            raiz = approach[iteration-1];
        }else{
            cout <<"\nRAIZ : "<<raiz<<endl;
        }
        char decide;
            cout << " |b|BISECCION\n |f| FALSAS\n |n|NEWTON :"<<endl;
            cin >> decide;
                iteration=0;
                rank_left=iz;
                rank_right=der;
                fin = false;
                for(int i=0;i<100;i++){approach[i]=0;}
                if(decide=='b'){calcular('b');}
                if(decide=='f'){calcular('f');}
                if(decide=='n'){cal_newton(raiz,1);}
        return raiz;
    }
    return 0;
}

int Raices::trunca(float dato,int ceros){
    int retorna = dato*pow(10,ceros); //10 decimales
    return retorna;
}

double Raices::falsa(double rank_left,double rank_right){
    return rank_right - (f(rank_right)*(rank_right-rank_left))/(f(rank_right)-f(rank_left));
}

double Raices::biseccion(double rank_left,double rank_right){
    return (rank_left+rank_right)/2;
}

double Raices::f(double unknown){
    double sum=0;int auxi=exponencial;
    for(int i=0;i<auxi;i++){
        sum+=coefficient[i]*pow(unknown,exponencial);
        exponencial--;
    }
    exponencial=auxi;
    return sum + constant;
}
    //methods Determinant------------

bool Determinant::divergir(float matriz[][6],int order){
        int x=0,y=0,bad=0;
        float directriz[order];
        for(x=0;x<order;x++){
            for(y=0;y<order;y++){
                if(abs(matriz[x][x])<abs(matriz[y][x])){
                    return false;
                }
            }
        }
        return true;
}

float Determinant::determinant(float matriz[][6],int order){
        float det = 0;int y=0;
        if(order==1){
            det = matriz[0][0];
        }else{
            for(y=0;y<order;y++){
                det += matriz[0][y]*cofactors(matriz,order,0,y);
            }
        }
        return det;
}

float Determinant::cofactors(float matriz[][6],int order,int row,int column){
        float submatriz[6][6];
        int order2 = order-1,x=0,y=0,f=0,c=0;
        for(x=0;x<order;x++){
            for(y=0;y<order;y++){
                if(x!=row && y!=column){
                    submatriz[f][c] = matriz[x][y];
                    c++;
                    if(c>=order2){
                        c=0;f++;
                    }
                }
            }
        }
        return pow(-1,row+column)*determinant(submatriz,order2);
}

void Determinant::get_data(){
        int x=0,y=0,v=0; char opc;
        do{
            cin.clear();
            fflush(stdin);
            cout << " [TAMANIO DE LA MATRIZ] : ";
            v=scanf("%d",&order);
        }while(v==0||order>6||order<2);
        do{
                opc='0';
            for(x=0;x<order;x++){
                for(y=0;y<order;y++){
                    do{
                        cin.clear();
                        fflush(stdin);
                        cout << " COEFICIENTES ["<<x<<"]["<<y<<"] : ";
                        v=scanf("%f",&matriz[x][y]);
                    }while(v==0);
                }
            }
            if(!divergir(matriz,order)){
                cin.clear();
                fflush(stdin);
                cout << " LA DIAGONAL PRINCIPAL NO TIENE LOS VALORES PREDOMINANTES !!! PUEDE DIVERGIR !!!"<<endl;
                cout << "           |r| PARA REACER LA MATRIZ OTRA TECLA PARA CONTINUAR"<<endl;
                cin >> opc;
            }
        }while(opc=='r');
            cout << " [TERMINOS INDEPENDIENTES] " << endl;
            for(x=0;x<order;x++){
                do{
                    cin.clear();
                    fflush(stdin);
                    cout << " ELEMENTOS X_"<<x+1<<" : ";
                    v=scanf("%f",&independent[x]);
                }while(v==0);
            }
}

    // methods Gauss ------------------
void Gauss::eliminacion(float matriz[][6],float independent[],int order){
        int j=0,i=0,k=0;
        float solucion[order],suma=0,aux=0;
        for(j=0;j<order;j++){
            solucion[j]=0;
        }
        // Triangulo -----
        //A21 = A21-(A21/A11)*A11
        //A22 = A22-(A21/A11)*A12
        //A23 = A23-(A21/A11)*A13
        //A31 = A31-(A31/A11)*A11
        //A32 = A32-(A31/A11)*A12
        for(k=0;k<order-1;k++){ //directriz (toda la columna abajo)
            for(i=k+1;i<order;i++){
                aux=matriz[i][k]/matriz[k][k];
                for(j=k;j<order;j++){
                    matriz[i][j]=matriz[i][j]-aux*matriz[k][j];
                }
                independent[i] -= aux*independent[k];
                print_matriz(matriz,order,independent);
            }
        }
        print_matriz(matriz,order,independent);

        // sustitucion hacia atras ---
        order--; // 2
        solucion[order] = independent[order]/matriz[order][order];
        for(i=order-1;i>=0;i--){
            suma=0;
            for(j=i+1;j<=order;j++){
                suma+=matriz[i][j]*solucion[j];
            }
            solucion[i] = (independent[i]-suma)/matriz[i][i];
        }

        cout << "  SOLUCION DE INCOGNITAS "<<endl;
        for(j=0;j<=order;j++){
            cout << "  x_"<<j<<": "<<solucion[j]<<"     ";
        }
}

void Gauss::jordan(float matriz[][6],float independent[],int order){
        float inversa[order][order],pivote=0,aux=0;
        int i=0,j=0,k=0;
        //inicializar inversa-----
        for(i=0;i<order;i++){
            for(k=0;k<order;k++){
                inversa[i][j]=0;
                if(i==j){inversa[i][j]=1;}
            }
            inversa[i][i]=1;
        }
        //gauss jordan process ---
        for(i=0;i<order;i++){
            pivote = matriz[i][i];
            independent[i] = independent[i]/pivote;
            for(k=0;k<order;k++){ //dividir la fila del pivote que volvimos 1
                matriz[i][k]=matriz[i][k]/pivote;
                inversa[i][k]=inversa[i][k]/pivote;
            }
            for(j=0;j<order;j++){ //volver cero lo que hay bajo y arriba del pivote
                if(i!=j){
                    aux=matriz[j][i];
                    independent[j] = independent[j]-aux*independent[i];
                    for(k=0;k<order;k++){
                        // 3        =   3        - elemento a cero * elemento de la fila pivote
                        matriz[j][k]=matriz[j][k]-aux*matriz[i][k];
                        inversa[j][k]=inversa[j][k]-aux*inversa[i][k];
                    }
                }
            }
        }
        print_matriz(matriz,order,independent);
        cout<<"\n\n : [INVERSA] : "<<endl;
        for(i=0;i<order;i++){
            for(k=0;k<order;k++){
                 if(matriz[i][k]>=0){cout <<"    "<<inversa[i][k];}else{cout<<"   "<<inversa[i][k];}
            }
            cout <<"\n";
        }
}

void Gauss::seidel(float matriz[][6],float independent[],int order){
        float incognitas[order],suma_fila=0,error_actual=0,aproximacion[order];
        int directriz=0,suma=0,iteracion=0,bucle=0;
        for(iteracion=0;iteracion<order;iteracion++){
            incognitas[iteracion]=0;
            aproximacion[iteracion]=0;
        }
        printf("\n         x1              x2                x3                ERROR :\n");
        do{ //cuantas iteraciones se haran
            for(directriz=0;directriz<order;directriz++){
                suma_fila=0;
                for(suma=0;suma<order;suma++){
                    if(suma!=directriz){
                        suma_fila = suma_fila+(-1*matriz[directriz][suma]*incognitas[suma]);
                    }
                }
                incognitas[directriz] = (independent[directriz] + suma_fila)/matriz[directriz][directriz];
            }
            for(iteracion=0;iteracion<order;iteracion++){
                error_actual = ((incognitas[iteracion]-aproximacion[iteracion])/incognitas[iteracion])*100;
                if(incognitas[iteracion]>=0){printf("|    %.6f    |",incognitas[iteracion]);}else{printf("|   %.6f    |",incognitas[iteracion]);}
            }printf("   %.6f\n",error_actual);
            for(iteracion=0;iteracion<order;iteracion++){
                if(aproximacion[iteracion]!=incognitas[iteracion]){
                    aproximacion[iteracion] = incognitas[iteracion];
                }else{
                    bucle++;
                }
            }
            if(bucle==order){break;}
            if(bucle>100||incognitas[0]-incognitas[3]>1000){cout << " !!!A DIVERGIDO!!!"<<endl;break;}
        }while(true);
}

// functions -------------

int main(){
    char opcion;
    cin.clear();
    fflush(stdin);
    cout <<" |L| SISTEMA DE ECUACIONES LINEALES\n |C| RAIZ DE LA ECUACION\n |s| SALIR"<<endl;
    cin >> opcion;
    switch(opcion){
        case 'l':
            lineales();
            break;
        case 'c':
            cuadraticas();
            break;
        case 's':
            return 0;
    }
    main();
}

void cuadraticas(){
    char method;
    Raices raiz;
    cin.clear();
    fflush(stdin);
    cout << " METODO\n |b|BISECCION\n |f|FALSA\n |n| NEWTON "<<endl; cin >>method;
    if(method=='b'||method=='f'){
        raiz.get_rank_left();
        raiz.get_rank_right();
        raiz.get_tolerance();
        raiz.get_polinomio();
    }
    if(method=='b'||method=='f'){
        raiz.calcular(method);
    }
    if(method=='n'){
        raiz.cal_newton(0,0);
    }
}

void lineales(){
    char opcion;
    float reemplazo[6][6],reemplazo_indep[6]; //matriz para reiniciar la matriz tras cada metodo
    Determinant deter;
    deter.get_data();
    for(int i=0;i<6;i++){
        for(int x=0;x<6;x++){
            reemplazo[i][x] = deter.matriz[i][x];
        }
        reemplazo_indep[i]=deter.independent[i];
    }
    do{
    print_matriz(deter.matriz,deter.order,deter.independent);
    float resultado = deter.determinant(deter.matriz,deter.order);
    if(resultado!=0){
        cout << "\n [DETERMINANTE] : "<< resultado << endl;
        cout << "\n |g| GAUSS |j| GAUSS-JORDAN |s| SEIDEL : " << endl;
        cin >> opcion;
        Gauss gauss;
        switch(opcion){
            case 'g':
                gauss.eliminacion(deter.matriz,deter.independent,deter.order);
                break;
            case 'j':
                gauss.jordan(deter.matriz,deter.independent,deter.order);
                break;
            case 's':
                gauss.seidel(deter.matriz,deter.independent,deter.order);
                break;
        }
    }else{
        cout << "\n DETERMINANTE ES INDETERMINADO O NO TIENE SOLUCION"<<endl;
    }
    stop;
    cls;
    cout << "\n  |c| PARA REALIZAR OTRO CALCULO\n  |n| PARA NUEVO SISTEMA DE ECUACIONES\n  |OTRO PARA SALIR|  "<<endl;
    cin >> opcion;
    for(int i=0;i<6;i++){
        for(int x=0;x<6;x++){
            deter.matriz[i][x]=reemplazo[i][x];
        }
        deter.independent[i]=reemplazo_indep[i];
    }
    }while(opcion=='c');
    if(opcion=='n'){
        main();
    }
}

void print_matriz(float matriz[][6],int order,float independent[]){
    int x=0,y=0;
     cout << "\n\n\n : [MATRIZ] : "<<endl;
    for(x=0;x<order;x++){
        for(y=0;y<order;y++){
            cout <<"\t"<<matriz[x][y];
        }
        cout << "|\t"<<independent[x]<<"\n";
    }
}












