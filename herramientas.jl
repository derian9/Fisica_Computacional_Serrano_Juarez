
__precompile__() # Este comando es para que julia precompile el paquete

module herramientas

export método_newton
function método_newton(f,df,s)
    x=s;
    for i in 1:20
    
    x=x-f(x)/df(x)
        
end
return x
end

export método_newton2
function método_newton2(f,df,inicial)
      
    list=zeros(Complex64,20);
    x=inicial;
for i in 1:20
        x=x-f(x)/df(x)
    list[i]=x; # Me crea una lista de datos que ha iterado
    end
    return x
end


export método_newtonlinspace
function método_newtonlinspace(f,df,a,b)
    x=a:b #ESte es nuestro linspace
    s=b-a #Es el numero de valores que va a realizar
    list=zeros(s)
for i in 1:s
    list[i]=método_newton2(f,df,x[i])
    end
    return list
end

export método_newton3
function método_newton3(f,df,a,b)
  
list=método_newtonlinspace(f,df,a,b) #Que os haga una lista usando el mètodo anterior 
   t=[]
    
    for i in 1:length(list)
    iteraciones=0
          for j in 1:length(t)
              if list[i]>t[j]+e
                    iteraciones=iteraciones+1
                end
                if list[i]<iteraciones[j]-e
                    iteraciones=iteraciones+1
                end
            end
            if iteraciones==length(t)
                push!(t,list[i]) #Le pedimos que nos que todo lo repetido
            end
        end
        return t
    end


export metodo_bieccion
function biseccion(f,a,b)
    p=(a+b)/2
for i in 1:20
        if  f(a)*f(p) > 0
            a = p
        end
        if f(b)*f(p) > 0
            b = p
        end
        p = (a+b)/2
    end
    return p,f(p)
end

export integral_rectangulo
function int_rectangulo(f,a,b,n) #Llamamos la funcion int con los armumentos deseados (f la función; a,b es el intervalo , n es la presicion con la que se va a calcular la integral.)
   
    s=0    #S será la variable que irá cargando la suma de los rectangulos. innicialmente en cero.
    dx=(b-a)/n #Aqui definimos la base de nuestros rectangulos.
    
    for i in 0:n
    
        s = s + dx*f((a+(2i+1)*dx)/2) #Damos la aproximacion de los trapecios.
    end
    return s/100 #Pedimos que regrese s, es decir guarde el valor de s y vuelva a realizar el procedimiento.
    
end

export integral_trapecio
function int_trapecio(f,a,b,n) #Llamamos la funcion int con los armumentos deseados (f la función; a,b es el intervalo , n es la presicion con la que se va a calcular la integral.)
   
    s=0    #S será la variable que irá cargando la suma de los rectangulos. innicialmente en cero.
    dx=(b-a)/n #Aqui definimos la base de nuestros rectangulos.
    
    for i in 0:n
    
        s = s + dx*(f(a+(i*dx))+f(a+((i+1)*dx)))/2  #Damos la aproximacion de los trapecios
    end
    return s/100
    
end

export intergral_simpson
function int_simpson(f,a,b,n)
  s=0
    dx= (b-a)/n;
    adx=a;
    bdx=adx+n;
    
    for i in 1:dx
        s+=(n/6)*(f(adx)+4*f((adx+bdx)/2)+f(bdx));
        adx = bdx
        bdx = bdx+n
    end
    
    return s
    
end


export derivada_númerica
function derivada_númerica(f,x,h) #Nombramos nuestra función "derivada númerica" con sus respectivas variables.
    
       df=(f(x+h)-f(x))/h #Introducimos la noción de derivada (la fórmula anteriormente vista).
    
end

export derivada_númerica2
function derivada_númerica2(f,x)
list=zeros(100) #Se crea una lista de ceros donde se guardará el valor de la derivada
n=1
           for i in 1:100 
              df=(f(x+(1/n))-f(x))*n #Aplicamos la derivada
              list[i]=df; #Asignamos el valor de la derivada a un elemento de la lista
              n+=1 # se aumenta el valor de n en cada iteración,aqui h=1/n
           end
return list 
end


export derivada_simetrica
function derivada_simetrica(f,x,h) #Creamos una función que toma como entradas la función , el punto a evaluar,h
              df=(f(x+h)-f(x-h))/2h #Aplicamos la fórmula
        
return df
end



export derivada_simetrica2
function derivada_simetrica2(f,x) #Realizamos lo mismo que en el ejercicio 2.
list=zeros(100)
n=1
           for i in 1:100
              df=(f(x+(1/n))-f(x-(1/n)))*n/2 
              list[i]=df; #Asignamos el valor de la derivada a un valor de la lista
              n+=1 #Se aumenta el valor de n en cada iteración,aqui h=1/n
           end
return list
end


export método_euler
function método_euler(f,x0,t0,tf,h)
    x=x0
    t=t0
    listx = [] #Creamos listas para colo car los valores y poder graficar posteriormente.
    listt = []
    push!(listt,t)
    push!(listx,x)
    n = round((tf-t0)/h) # Aqui aplicamos la fórmula de Euler.
    for i in 1:n-1
        x += h*f(x,t)
        t += h
        push!(listt,t)
        push!(listx,x)
    end
    return listt, listx
end


export método_euler2
function método_euler2(g,x0,u0,listt)
     n=length(listt)
     listx=zeros(n)             # Creamos una lista para las x|k del método de Euler.
     listu=zeros(n)             # Creamos una lista para las u|k.
     listx[1]=x0                # Definimos las condiciones iniciales.
     listu[1]=u0                
     h=(listt[n]-listt[1])/n   #Aplicamos la fórmula de Euler.
        for i in 1:n-1
        listx[i+1]=listx[i]+h*listu[i]                     
        listu[i+1]=listu[i]+h*g(listx[i],listu[i],listt[i])
       end
     return listx  #Pedimos que nos devuelva x para poder gráficsr.
end


export método_euler3
function método_euler3(h,x0,listt) #Creamos una función para aproximar, de la misma forma que antes.
    n = length(listt)            
    listx = zeros(n)             
    listx[1] = x0               
    m = (listt[n]-listt[1])/n    
    for i in 1:n-1
        listx[i+1] = listx[i] + m*h(listx[i],listt[i])   #Aplicamos la fórmula del método de Euler.
    end
    return listx                  
end


export met_newton
function met_newton(g,x0,h) #Definamos la función de aproximacion de Newton como g(x)=x-f(x,t_1)h.
    dg(z) = derivada_simetrica(g,z,h) #Aplicando para la función g, x0 y h.
    x = x0;                            
    for i in 1:10                     
        x = x-g(x)/dg(x) #Usamos la fórmula de recurrencia del método deNewton. 
    end
    return x                          
end


export punto_fijo
function punto_fijo(f,x0)
    x = x0         
    for i in 1:20
        x = f(x) # Iteramos para aproximar el punto fijo.
    end
    return x      
end



export EulerImplicito
function EulerImplicito(f,x0,listt,metodo)
    n = length(listt)   # N° elementos de la lista listt.
    listx = zeros(n)                   
    listx[1] = x0       #Colocamos la condición inicial x=x0.
    h = listt[2]-listt[1]                
    for k in 1:n-1
        xk = listx[k]                    
        t = listt[k+1]                  
        if metodo == "newton"
            g(z) = z - xk - h*f(z,t)     #g es la función a la que queremos encontrar su raíz x_(k+1).
            listx[k+1] = met_newton(g,xk,h)  
        elseif metodo == "fijo"
            G(z) = xk + h*f(z,t)
            listx[k+1] = putno_fijo(G,xk)     
       
        end
    end
    return listx                         
end


export Punto_Medio
function Punto_Medio(f,listt,inicial)
    h=listt[2]-listt[1] #Esta es a distancia entre las primeras dos entradas.
    listx=zeros(length(listt)) #Misma longitud de las listas.
    listx[1]=inicial #Apliquemos las iteraciones dada la condición de recurrencia.
        for i in 1:length(listt)-1
           listx[i+1]= listx[i]+h*f(listx[i]+(h/2)*f(listx[i],listt[i]),listt[i]+h/2)
        end 
    return listx #Pedimos que regrese los resultados.
end



export rungeKutta4
function rk4(f,list,x0) #De igual forma con Runge-Kutta de orden 4.
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        k1 = f(x,t)
        k2 = f(x+(h/2)*k1,t+(h/2))
        k3 = f(x+(h/2)*k2, t+(h/2))
        k4 = f(x+h*k3, t+h);
        x = x+(h/6)*(k1+2*k2+2*k3+k4)
        push!(listx,x) 
     end
     return listx
end




end


import herramientas.jl
