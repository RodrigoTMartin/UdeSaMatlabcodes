clear;
clc;
close all;

plot_options.model_names={'TNT_NK_2_agasto','TNT_NK_2_baseline1gasto'};% Nombre de los .mat creados al final del correr Dynare correspondientes a los distintos modelos
shocks        ={'Rw' };%nombre de variables shockeadas, con el mismo nombre que el declarado en Dynare (este programa asume que si la vraiable es, por ejemplo, 'a', el shocks es llamado 'u_a' en el mod-file)
shocks_latex  ={'R^W'};  %nombre que aparecera en el titulo de cada IRF para cada variable ex?gena. Se permite escribir en formato LaTex

v_sel                   = {'c' ,'cW' ,'h','d' ,'cC','hW','tao','cT','cN'}; %nombre de las variables cuyo irf desean mostrarse, con el mismo nombre que el declarado en Dynare
v_name                  = {'c','c^W','h','d','c^C','h^W','tao','c^T','c^N'};%nombre que aparecer? en el titulo de cada IRF para cada variable ex?gena. Se permite escribir en formato LaTex.
plot_options.v_do_cusum = [    0,    0,    0,    0,      0,     0,    0,  0];%Vector (de tama?o igual al numero de variables) con un 1 en la posicion de la variable en la cual se desea desplegar la IRF acumulada (util cuando la variable es una tasa de crecimiento, por ejemplo se usa para graficar el tipo de cambio S en lugar de la tasa de depreciacion ppiS, y M en lugar del crecimiento del dinero d_M). 
plot_options.v_do_cusum = [plot_options.v_do_cusum 0]; % Se agrega un elemento porque tambien se grafica la IRF de la variable exogena. 
plot_options.v_adj      = 100*ones(1,size(v_sel,2)+1); % Vector (de tama?o igual al numero de variables +1) con factores de ajuste para cada variable. En este caso, se multiplican todas por 100 para tener las respuestas en puntos porcentuales. Se agrega un elemento (+1) porque tambien se grafica la IRF de la variable exogena.
plot_options.v_div_ss   = zeros(1,size(v_sel,2)+1); % Vector (de tama?o igual al numero de variables +1) con unos en caso de quere dividr la respuesta por el valor de estado estacionario de cada varaible (util si las aproximacion es lineal en lugar de log-lineal). En este caso son todos ceros pues esto no lo usaremos. Se agrega un elemento (+1) porque tambien se grafica la IRF de la variable exogena.
plot_options.horizon    = 20; % numero que indica el horizonte al mostrar de la IRF, que debe ser menor o igual al numero en la opcion irf de stcoh_simul en Dynare.
plot_options.n_col      = 4; %numero de columnas en el sub-plot.
plot_options.n_row      = 4; %numero de filas en el sub-plot. El  numero de filas multiplicado por el de columnas debe ser mayor o igual al numero de variables que quieren ser graficadas +1 (recordando que el programa siempre grafica la respuesta de la variable que esta siendo shockeada).
plot_options.grid       = 1; % 1 si se incluye una grilla de puntos en el grafico
plot_options.latex      = 1; % 1 si los titulos son escritos con formato Latex

plot_options.saving=0; % opcion que permite guardar cada uno de los grafico en el current folder. Si es 0 no guarda nada, si se 1 guarda en formato eps, si es 2 en formato pdf, si es 3 en formato jpeg. El nombre del archivo creado esta identificado por el nombre del shock.
plot_options.marks={'-b','--r','-.k',':m'}; % opcion que permite guardar cada uno de los grafico en el current folder. Si es 0 no guarda nada, si se 1 guarda en formato eps, si es 2 en formato pdf, si es 3 en formato jpeg. El nombre del archivo creado esta identificado por el nombre del shock.

% De aqui en mas no tocar nada.


for j=1:size(shocks,2)
    plot_options.u_sel=shocks(j);
    if size(plot_options.model_names,2)>1
        for i=2:size(plot_options.model_names,2)
            plot_options.u_sel=[plot_options.u_sel shocks(j)];
        end
    end
    plot_options.v_sel=[v_sel  shocks(j)];
    plot_options.u_name=shocks_latex(j);
    plot_options.v_name=[ v_name shocks_latex(j)];
    plots_for_dynare(plot_options)
    if plot_options.saving==1
        saveas(gcf,['irf_' shocks{j} '.eps'],'psc2');
    elseif plot_options.saving==2
        saveas(gcf,['irf_' shocks{j} '.pdf'],'pdf');
    elseif plot_options.saving==3
        saveas(gcf,['irf_' shocks{j} '.jpg'],'jpg');
    end
end
