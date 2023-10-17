
#Definição das condições de contorno:
#Condição C1, próximo ao estado trivial de equilíbrio da linearização:
C1 = [0.1 0.1 0 0];
t1 = 0:0.001:100;
#Condição C2, afastada do equilíbrio trivial:
C2 = [(pi()/3) (-pi()/4) 0 0];
t2 = 0:0.001:30;

#Funções de representação em Espaço de Estados
#Não-linearizada:
function dx1 = f1(t,x)
    wp = 1; %omega_p
    lda = 1.9; %lambda
    dx1(1) = x(3);
    dx1(2) = x(4);
    dx1(3) = -(3*wp^2*((5+4*lda)*sin(x(1)) + 3*sin(x(1) - 2*x(2))))/(15+8*lda-9*cos(2*(x(1)-x(2))) - (9*sin(2*(x(1)-x(2))*x(3)^2))/(15+8*lda-9*cos(2*(x(1)-x(2)))) - (6*sin(x(1)-x(2))*x(4)^2)/(lda*(4*(3+4*lda)-9*(cos(x(1)-x(2)))^2)));
    dx1(4) = -(3*wp^2*(-3*(2+lda)*sin(2*x(1)-x(2))+(6+lda)*sin(x(2))))/(15+8*lda-9*cos(2*(x(1)-x(2)))) + (6*lda*(3+lda)*sin(x(1)-x(2))*x(3)^2)/(4*(3+lda)-9*(cos(x(1)-x(2)))^2) + (9*sin(2*(x(1)-x(2)))*x(4)^2)/(15+8*lda-9*cos(2*(x(1)-x(2))));
end

#Linearizada:
function dx2 = f2(t,x)
    wp = 1; %omega_p
    lda = 1.9; %lambda
    dx2(1) = x(3);
    dx2(2) = x(4);
    dx2(3) = -(3*wp^2*(2*(2+lda)*x(1)-3*x(2)))/(3+4*lda);
    dx2(4) = (3*lda*(3*(2+lda)*x(1) - 2*(3+lda)*x(2)))/(3+4*lda);
end

#Integração Numérica
#Método M1 = Dormand-Prince (C1) :
[t1_lin1, x1_lin1] = ode45(@f2, t1, C1); %Linearizada
[t1_nonlin1, x1_nonlin1] = ode45(@f1, t1,C1); %Não-Linearizada
#Método M2 = Método Bogacki-Shampine (C1):
[t1_lin2, x1_lin2] = ode23(@f2, t1, C1); %Linearizada
[t1_nonlin2, x1_nonlin2] = ode23(@f1, t1, C1); %Não-Linearizada

#Método M1 = Dormand-Prince (C2) :
[t2_lin1,x2_lin1] = ode45(@f2, t2, C2); %Linearizada
[t2_nonlin1,x2_nonlin1] = ode45(@f1, t2,C2); %Não-Linearizada
#Método M2 = Método Bogacki-Shampine (C2):
[t2_lin2, x2_lin2] = ode23(@f2, t2, C2); %Linearizada
[t2_nonlin2, x2_nonlin2] = ode23(@f1, t2, C2); %Não-Linearizada


#Análise Gráfica sobre a condição C1
#Usando o método M1:
figure(1)
plot(x1_lin1(:,3), x1_lin1(:,1), "r");
hold on;
plot(x1_nonlin1(:,3), x1_nonlin1(:,1), "b");
legend('Linearizado','Não-linearizado');
xlabel('d\theta_1/dt')
ylabel('\theta_1')
title('Comparação entre \theta_1 em função de d\theta_1/dt na integração das equações linearizadas e não-linearizadas')

figure(2)
plot(x1_lin1(:,4), x1_lin1(:,2), "m");
hold on;
plot(x1_nonlin1(:,4), x1_nonlin1(:,2), "g");
legend('Linearizado','Não-linearizado');
xlabel('d\theta_2/dt')
ylabel('\theta_2')
title('Comparação entre \theta_2 em função de d\theta_2/dt na integração das equações linearizadas e não-linearizadas')

#Usando o método M2:
figure(3)
plot(x1_lin2(:,3), x1_lin2(:,1), "r");
hold on;
plot(x1_nonlin2(:,3), x1_nonlin2(:,1), "b");
legend('Linearizado','Não-linearizado');
xlabel('d\theta_1/dt')
ylabel('\theta_1')
title('Comparação entre \theta_1 em função de d\theta_1/dt na integração das equações linearizadas e não-linearizadas')

figure(4)
plot(x1_lin2(:,4), x1_lin2(:,2), "m");
hold on;
plot(x1_nonlin2(:,4), x1_nonlin2(:,2), "g");
legend('Linearizado','Não-linearizado');
xlabel('d\theta_2/dt')
ylabel('\theta_2')
title('Comparação entre \theta_2 em função de d\theta_2/dt na integração das equações linearizadas e não-linearizadas')

#Análise Gráfica sobre a condição C2
#Usando o método M1:
figure(5)
plot(x2_nonlin1(:,3), x2_nonlin1(:,1), "b");
legend('Não-linearizado');
xlabel('d\theta_1/dt')
ylabel('\theta_1')
title('\theta_1 em função de d\theta_1/dt na integração do modelo não-linearizado')

figure(6)
plot(x2_nonlin1(:,4), x2_nonlin1(:,2), "g");
legend('Não-linearizado');
xlabel('d\theta_2/dt')
ylabel('\theta_2')
title('\theta_2 em função de d\theta_2/dt na integração do modelo não-linearizado')

#Usando o método M2:
figure(7)
plot(x2_nonlin2(:,3), x2_nonlin2(:,1), "b");
legend('Não-linearizado');
xlabel('d\theta_1/dt')
ylabel('\theta_1')
title('\theta_1 em função de d\theta_1/dt na integração do modelo não-linearizado')

figure(8)
plot(x2_nonlin2(:,4), x2_nonlin2(:,2), "g");
legend('Linearizado','Não-linearizado');
xlabel('d\theta_2/dt')
ylabel('\theta_2')
title('\theta_2 em função de d\theta_2/dt na integração do modelo não-linearizado')

#Parãmetros para cálculo da Energia Mecânica:
lda = 1.9; %adim.
wp = 1; %rad/s
g = 9.81; %m/s^s
mu = 0.2; %kg/m
l_1 = g/wp^2;
l_2 = l_1/lda;
m1 = mu*l_1;
m2 = mu*l_2;

#Energia Mecânica para M1:
T = (1/6)*(l_1^2)*(m1 + 3*m2)*(x2_nonlin1(:,3).^2) + (1/2)*cos(x2_nonlin1(:,1)-x2_nonlin1(:,2))*l_1*l_2*m2*(x2_nonlin1(:,3))*(x2_nonlin1(:,4)) + (1/2)*(l_2^2)*m2*(x2_nonlin1(:,4).^2);
V = -(1/2)*m2*g*l_2*cos(x2_nonlin1(:,2)) - (1/2)*(m1+2*m2)*g*l_1*cos(x2_nonlin1(:,1));
M = T + V;
figure(9)
plot(t2_nonlin1, T, 'c');
hold on;
plot(t2_nonlin1, V, 'g');
hold on;
plot(t2_nonlin1, M, 'b');
legend('Energia Cinética','Energia Potencial', 'Energia Mecânica');
xlabel('Tempo (s)');
ylabel('Energia (J)');
title('Variação da Energia Mecânica com o tempo no método de integração M1');


