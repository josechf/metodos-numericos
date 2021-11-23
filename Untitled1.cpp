
#include <bits/stdc++.h>

#define stop system("pause");
#define cls system("cls");
#define buff fflush(stdin);

//utilizar el rango como en integral en newton, lagrange, minimos, para que sea mas optimo

// nota las funciones get deben ser set, pero me da pereza cambiarlas
// nota los atributos deberian ser private, pero me da pereza cambiarlos
// creo que hay vainas sin validar, pero yo que se

using  namespace std;
void print_matriz(float[][6],int,float[]);
void lineales();
void cuadraticas();
void interpolacion();
void diferencial();
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

class Interpolation{
    public:
        int amount,x,y,range[2];
        float **date_table,unknow,**sub_matriz;
        int n;
    public:
        Interpolation();
        void exchange();
        void print_matriz(float**,int);
        void set_dates();
        void set_range();
        void set_unknow();
        void newton();
        void regresion();
        void cuadro_suma(float[]);
        void lagrange();
        void integral();
        float simpson13();//simpson 1/3
        float simpson38(); // simpson 3/8
        float trapesoidal();
        void calc_rango_h(int);
        void set_rango_n(int);   //pido valor de n para calcular h
        bool secuencia();    //esto es bool para retornar false en caso que el rango h este mal
};

class Diferenciales{
    public:
        float *y,*k,x,h,Y,*coeficiente,fin;
        int maxim,elem,*opcion,c;
    public:
        float euler(float,float);
        void euler_mejorado();
        float f(float,float);
        void set_data();
        void runge_kutta();
        void set_polinomio();
        float ec(int,float,float);
};

    // DIFERENCIAL

void Diferenciales::runge_kutta(){
    int a=0,i=0; float sum_k=0;
    set_data();
    cout<<"\tn\tx\ty";
    for(i=1;i<=maxim;i++){
        cout<<"\t\tk"<<i;
        k[i]=0;
    }
    cout<<"\t\tYn+1"<<endl;
    for(i=0;i<maxim;i++){
        cout<<"\t"<<i<<"\t"<<x<<"\t"<<Y<<"\t";
        k[0] = f(x,Y);
        sum_k = k[0];
        cout<<"\t"<<k[0]<<"\t";
        k[1] =  f(x+(0.5*h),Y+(0.5*k[0]*h));
        cout<<"\t"<<k[1]<<"\t";
        k[2] =  f(x+(0.5*h),Y+(0.5*k[1]*h));
        cout<<"\t"<<k[2]<<"\t";
        k[3] =  f(x+h,Y+(k[2]*h));
        cout<<"\t"<<k[3]<<"\t";
        Y += (h*(k[0]+2*k[1]+2*k[2]+k[3]))/6;
        cout<<"\t"<<Y<<endl;
        x+=h;
    }
}

void Diferenciales::set_data(){
    int i=0;
    cout<<" VALOR INICIAL DE X"<<endl;cin>>x;
    cout<<" VALOR INICIAL DE Y"<<endl;scanf("%f",&Y);
    cout<<" NUMERO DE ITERACIONES"<<endl;cin>>maxim;
    cout<<" DISTANCIA ENTRE DATOS"<<endl;cin>>h;
    set_polinomio();
    y = (float*) malloc(maxim);
    y[0]=Y;
    k = (float*)malloc(maxim-1);
    k[0] = 0;
    for(i=1;i<maxim;i++){
        y[i] = 0;
        k[i] = 0;
    }
}

void Diferenciales::euler_mejorado(){
    int i=0;
    set_data();
    cout<<"\tN\t"<<"Xn\t"<<"Yn\t"<<endl;
    cout<<"\t"<<0<<"\t"<<x<<"\t"<<y[0]<<endl;
    for(i=1;i<maxim;i++){
        y[i]=euler(x,y[i-1]);
        x+=h;
        y[i] = y[i-1]+h*(f(x-h,y[i-1])+f(x,y[i]))/2;
        cout<<"\t"<<i<<"\t"<<x<<"\t"<<y[i]<<endl;
    }
}

float Diferenciales::euler(float x, float y){
    return y + h*(f(x,y));
}
void Diferenciales::set_polinomio(){
    int i=0;
    cout <<" CANTIDAD DE ELEMENTOS : ";cin>>elem;
    coeficiente = (float*) malloc(elem);
    opcion = (int*) malloc(elem);
    for(i=0;i<elem;i++){
        cout << " INGRESE ELEMENTO : "<<endl;
        //cout<<" 1) A\t 2) BX\t 3) CY\t 4) X^#\t 5) Y^#\n 6) BX*CY\t 7) DX/EY\n\n";cin>>opcion[i];
        cout<<" 1) A\t 2) BX\t 3) CY\t 4) X^#\t 5) Y^#\n";cin>>opcion[i];
    }
    for(i=0;i<elem;i++){
        if(opcion[i]==7||opcion[i]==6){
            cout<<" COEFICIENTE X";cin>>coeficiente[i];
            i++;
            cout<<" COEFICIENTE Y";cin>>coeficiente[i];
        }else{
        cout<<" COEFICIENTE : ";cin>>coeficiente[i];
    }
    }
}

float Diferenciales::f(float x,float y){ //debo cambiarlo por una funcion
    //return x-y+5;
    //return x*sqrt(y);
    c=0;float sum=0;int m=0;
    for(c=0;c<elem;c++){
        sum += ec(m,x,y);m++;
    }
    return sum;
}
float Diferenciales::ec(int m,float x,float y){
        if(opcion[m]==1){
            return coeficiente[c];
        }
        if(opcion[m]==2){
            return coeficiente[c]*x;
        }
        if(opcion[m]==3){
            return coeficiente[c]*y;
        }
        if(opcion[m]==4){
            return pow(x,coeficiente[c]);
        }
        if(opcion[m]==5){
            return pow(y,coeficiente[c]);
        }
        if(opcion[m]==6){
            return (coeficiente[c]*x)*(coeficiente[c+1]*y);c++;
        }
        if(opcion[m]==7){
            return (coeficiente[c]*x)/(coeficiente[c+1]*y);c++;
        }
}

    // INTERPOLACION

Interpolation::Interpolation(){
    amount=x=y=range[0]=range[1]=unknow=n=0;
}

void Interpolation::set_unknow(){
    do{
        cout << " VALOR DE INCOGNITA : "<<endl;
        x=scanf("%f",&unknow);
    }while(x==0);
}
void Interpolation::set_range(){
    cout << "\n | RANGO | " <<endl;int p=0;
    do{
        cin.clear();buff;
        cout<<" inicial : ";p=scanf("%d",&range[0]);
    }while(p==0||range[0]<0||range[0]>amount-1);
    do{
        cin.clear();buff;
        cout<<" final : ";scanf("%d",&range[1]);
    }while(p==0||range[1]>amount-1||range[1]<1||range[1]==range[0]);
}
void Interpolation::set_dates(){
    int p=0;x=0;y=0;
    do{
        cin.clear();buff;
        cout << " CANTIDAD DE DATOS " << endl;
        p=scanf("%d",&amount);
    }while(p==0||amount<2);
    date_table = new float*[amount]; //filas
    do{
        date_table[x] = new float[2]; //columnas
        do{
            cin.clear();buff;
            printf(" X_%d : ",x);
            p=scanf("%f",&date_table[x][y]);y++;
        }while(p==0);
        do{
            cin.clear();buff;
            printf(" f(x_%d) : ",x);
            p=scanf("%f",&date_table[x][y]);x++;y=0;
        }while(p==0);
        cout<<endl;
    }while(x<amount);
}
void Interpolation::set_rango_n(int par_impar){
    int i=0,vueltas=0;float f=0;n=6;
    do{
        cin.clear();
        buff;
        vueltas=0;
        cout <<" SEGMENTOS DE AREA N : "<<endl;
        scanf("%d",&n);
        if(par_impar==2){ //n debe ser par
            i=n/2;
            f=(float)(n/2);
            if(n%2==1){
                vueltas=1;
            }
        }
        if(par_impar==3){ //n debe ser impar
            i=n/2;
            f=(float)(n/2);
            if(n%2==0){
                vueltas=1;
            }
        }
    }while(n<0||n>=amount||vueltas==1);
}

void Interpolation::exchange(){
    float matriz_aux[amount][2];
    for(x=0;x<amount;x++){
        matriz_aux[x][1] = date_table[x][0];
        matriz_aux[x][0] = date_table[x][1];
    }
    for(x=0;x<amount;x++){
        date_table[x][0] = matriz_aux[x][0];
        date_table[x][1] = matriz_aux[x][1];
    }
    print_matriz(date_table,amount);
}

void Interpolation::calc_rango_h(int par_impar){   //no se donde esta el problema pero no pasa de 6 cuando n=8
    set_rango_n(par_impar);
    float aux=date_table[range[0]][0];float h=0;y=1;
    h = (date_table[range[1]][0]-date_table[range[0]][0])/n;
    sub_matriz = new float*[n+1]; //filas
    sub_matriz[0] = new float[2]; //columna
    sub_matriz[0][0] = date_table[range[0]][0];
    sub_matriz[0][1] = date_table[range[0]][1];
    aux=aux+h;
    for(x=range[0];x<=range[1];x++){
        if(aux==date_table[x][0]){
            sub_matriz[y] = new float[2]; //columna
            sub_matriz[y][0] = date_table[x][0];
            sub_matriz[y][1] = date_table[x][1];
            y++;
            aux=aux+h;
        }
    }
    print_matriz(sub_matriz,n+1);
}

bool Interpolation::secuencia(){
    int aux_range=0;float h=0,aux=0;
    h=date_table[range[0]+1][0]-date_table[range[0]][0];
    aux=date_table[range[0]][0];
    for(x=range[0];x<=range[1];x++){
        if(aux!=date_table[x][0]){
            return false;
        }
        aux+=h;
    }
    return true;
}

void Interpolation::integral(){
    int opcion=0,aux=0;
    float resultado1=0,resultado2=0,resultado3=0;
    cout<<"\n 1.[TRAPESOIDAL]\n 2.[SIMPSON 1/3]\n 3.[SIMPSON 3/8]\n 4.[COMBINAR 1/3 3/8]";
    cin>>opcion;
    if(opcion==4){
        aux++;
        cout<<" INICIAR CON: 2.[1/3]\t3.[3/8]\n";
        cin>>opcion;
    }
    do{
        switch (opcion)
        {
        case 1:
            calc_rango_h(1);
            resultado1=trapesoidal();
            break;
        case 2:
            if(aux>0){set_range();}
            calc_rango_h(2);
            resultado2=simpson13();
            if(aux>0){opcion=3;aux++;}
            break;
        case 3:
            if(aux>0){set_range();}
            calc_rango_h(3);
            resultado3=simpson38();
            if(aux>0){opcion=2;aux++;}
            break;
        }
        if(aux>2){
            cout<<" SIMPSON COMBINADO : "<<resultado2+resultado3<<endl;
            aux=0;
        }
    }while(aux!=0);
}

float Interpolation::simpson13(){
    float suma_impar=0,suma_par=0,resultado=0,h = 0;
    h = (date_table[range[1]][0]-date_table[range[0]][0])/n;
    printf("\n\n SIMPSON 1/3 : ");
    for(x=1;x<n;x++){
        suma_impar += sub_matriz[x][1];x++;
    }
    for(x=2;x<n;x++){
        suma_par+=sub_matriz[x][1];x++;
    }
    resultado = h*((sub_matriz[0][1]+4*(suma_impar)+2*(suma_par)+sub_matriz[n][1])/3);
    cout<<resultado<<endl;
    return resultado;
}

float Interpolation::simpson38(){
    float suma_impar=0,suma_par=0,resultado=0,h = 0;
    printf("\n\n SIMPSON 3/8 : ");
    h = (date_table[range[1]][0]-date_table[range[0]][0])/n;
    for(x=1;x<n;x++){
        suma_par += 3*sub_matriz[x][1];
    }
    resultado = (h)*(sub_matriz[0][1]+suma_par+sub_matriz[n][1])/8;
    cout<<resultado<<"\n\n";
    return resultado;
}

float Interpolation::trapesoidal(){
    float suma_impar=0,suma_par=0,resultado=0,h = 0;
    printf("\n\n TRAPEZOIDE : ");
    //aqui suma_par la uso para sumar todo
    suma_par=0;
    for(x=1;x<n;x++){
        suma_par += sub_matriz[x][1];
    }
    resultado = (h/2)*(sub_matriz[0][1]+2*(suma_par)+sub_matriz[n][1]);
    cout<<resultado<<"\n\n"<<endl;
    return resultado;
}

void Interpolation::lagrange(){
    float G=0,renglon=0;
    int aux_rango=0;
    for(x=range[0];x<=range[1];x++){
        aux_rango++;
    }
    set_unknow();
    for(x=0;x<aux_rango;x++){
        renglon=1;
        for(y=0;y<aux_rango;y++){
            if(x!=y){
                renglon *= (unknow-date_table[y][0])/(date_table[x][0]-date_table[y][0]);
            }
        }
        G += date_table[x][1]*renglon;
    }
    cout << "\n\n G(x) = "<<G<< endl;
}

void Interpolation::newton(){
    // uso sub matriz para solo hacer newton con el rango que me dieron
    int h=0,aux_rango=0;
    float B=0,aux=0;
    for(x=range[0];x<=range[1];x++){
        aux_rango++;
    }
    float coefficent[aux_rango][aux_rango],sub_matriz[aux_rango][2];
    for(x=0;x<aux_rango;x++){
        for(y=0;y<aux_rango;y++){
            coefficent[x][y]=0;
        }
    }
    for(x=0;x<aux_rango;x++){
        for(y=0;y<2;y++){
            sub_matriz[x][y]=0;
        }
    }
    y=range[0]; // y es el rango para la tabla original, h es para ir de 0 a 1 en y, aux es la matriz original pero solo con los elementos del rango
    for(x=0;x<=aux_rango;x++){
        coefficent[x][0] = date_table[y][1];
        for(h=0;h<=2;h++){
            sub_matriz[x][h] = date_table[y][h];
        }
        if(y==range[1]){break;}
        y++;
    }
    for(y=1;y<aux_rango;y++){
        for(x=y;x<=range[1];x++){
            coefficent[x][y] = (coefficent[x][y-1]-coefficent[x-1][y-1])/(sub_matriz[x][0]-sub_matriz[x-y][0]);
        }
    }
    printf("\n TABLA \n"); //eliminar esta linea de codigo
    for(x=0;x<aux_rango;x++){
        for(y=0;y<aux_rango;y++){
            printf("\t%f ",coefficent[x][y]);
        }
        printf("\n");
    }
    for(x=0;x<aux_rango;x++){   //imprimir polinomio
        if(coefficent[x][x]>0){printf("+");}
        if(coefficent[x][x]!=0){
            printf("%f",coefficent[x][x]);
            for(y=0;y<x;y++){
                printf("(x-%f)",sub_matriz[y][0]);
            }
        }
    }
    set_unknow();
    for(x=0;x<aux_rango;x++){
        aux=1;
        for(y=0;y<x;y++){
            aux *= (unknow-sub_matriz[y][0]);
        }
        B += coefficent[x][x]*aux;
    }
    cout<<" RESULTADO : "<<B<<endl;
}

void Interpolation::cuadro_suma(float summation[]){
    float operacion[amount];
    for(x=0;x<amount;x++){
        summation[0] += date_table[x][0];
    }
    for(x=0;x<amount;x++){   //ln y
        operacion[x] = log(date_table[x][1]);
    }
    for(x=0;x<amount;x++){
        summation[1] += operacion[x];
    }
    for(x=0;x<amount;x++){   //x*ln y
        operacion[x] = date_table[x][0]*operacion[x];
    }
    for(x=0;x<amount;x++){
        summation[2] += operacion[x];
    }
    for(x=0;x<amount;x++){   //x^2
        operacion[x] = date_table[x][0]*date_table[x][0];
    }
    for(x=0;x<amount;x++){
        summation[3] += operacion[x];
    }
}

void Interpolation::regresion(){
    float summation[4],A0=0,A=0,A1=0,Y=0;
    int aux_rango=0;
    char option;
    for(x=range[0];x<=range[1];x++){
        aux_rango++;
    }
    for(x=0;x<4;x++){
        summation[x]=0;
    }
    cuadro_suma(summation);
    for(x=0;x<4;x++){
        printf("  %f  ",summation[x]);
    }
    A1 = (aux_rango*(summation[2])-summation[0]*summation[1])/(aux_rango*summation[3]-summation[0]*summation[0]);
    A0 = (summation[1]/aux_rango)-(A1*summation[0]/aux_rango);
    A = log(A0);
    A = exp(A0);
    cin.clear();
    fflush(stdin);
    cout<<" |r| REEMPLAZAR EN X\n";
    scanf("%c",&option);
    if(option!='r'){
        printf(" Y = %fe^%fx",Y,A,A1);
        return;
    }else{
        cin.clear();
        fflush(stdin);
        set_unknow();
        if(x==1){
            Y = A*exp(A1*unknow);
            printf(" %f = %fe^%f*%f",Y,A,A1,unknow);
        }
    }
}

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
    cout <<" |L| SISTEMA DE ECUACIONES LINEALES\n |C| RAIZ DE LA ECUACION\n |I| INTERPOLACION\n |d| DIFERENCIAL\n |s| SALIR"<<endl;
    cin >> opcion;
    switch(opcion){
        case 'l':
            lineales();
            break;
        case 'c':
            cuadraticas();
            break;
        case 'i':
            interpolacion();
            break;
        case 'd':
            diferencial();
            break;
        case 's':
            return 0;
    }
    main();
}


void diferencial(){
    char decide;
    Diferenciales dif;
    do{
        cout<<" |r| RUNGE KUTTA\n |e| EULER MEJORADO\n |s| SALIR"<<endl;
        cin>>decide;
        switch (decide)
        {
        case 'r':
            dif.runge_kutta();
            break;
        case 'e':
            dif.euler_mejorado();
            break;
        default:
            return;
        }
    }while(true);
}

void interpolacion(){
    int x=0,y=0,decide=0;bool revision;
    char opcion;
    Interpolation interpo;
    interpo.set_dates();
    while(true){
    interpo.print_matriz(interpo.date_table,interpo.amount);
    cin.clear();
    buff;
    cout << " |s| INTERCAMBIAR X POR Y : "<<endl;
    cin>>opcion;
    if(opcion=='s'){
        interpo.exchange();
    }
    cin.clear();
    buff;
    cout << " |n| NEWTON \n |r| REGRESION\n |l| LAGRANGE\n |i| INTEGRAL\n|s| SALIR"<<endl;
    cin >>opcion;
    cin.clear();
    fflush(stdin);
    switch(opcion){
        case 'n':
            interpo.set_range();
            interpo.newton();
        break;
        case 'r':
            interpo.set_range();
            interpo.regresion();
        break;
        case 'l':
            interpo.set_range();
            interpo.lagrange();
        break;
        case 'i':
            do{
                interpo.set_range();
                revision=interpo.secuencia();
                if(revision){
                    interpo.integral();decide=0;
                }else{
                    cout<<" LA SECUENCIA X ESTA MAL"<<endl;
                    do{
                        cout<<" 1.[NUEVO RANGO]\n 2.[NUEVA TABLA]"<<endl;
                        x=scanf("%d",&decide);
                    }while(x==0);
                    if(decide==2){
                        free(interpo.date_table);
                        interpo.set_dates();
                        decide=1;
                    }
                }
            }while(decide==1);
        break;
        case 's':
            return;
    }
    stop;
    }
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

void Interpolation::print_matriz(float** matriz,int amount){
    int x=0,y=0;
    for(x=0;x<amount;x++){
        cout <<x<<":   ";
        for(y=0;y<=1;y++){
            cout <<"\t"<< matriz[x][y];
        }
        printf("\n");
    }
}
