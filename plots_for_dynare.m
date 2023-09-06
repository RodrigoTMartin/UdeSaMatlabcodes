function plots_for_dynare(plot_options)

% This function creates a graph that combines the impulse responses of
% different models and/or different shocks for teh same model, generates by
% Dynare. If the responses come from different models, it is assumed that
% all the models have the same names for variables and shocks in Dynare.
% The input plot_options is a structure that has to include the following:
%
%   plot_options.model_names: A set of strings containing the mat-file names
%       (without includeing the .mat extension) were the oucomes from dynare 
%       have been stored. The function will create a line in the graph for 
%       each string in model_names.
%   plot_options.marks: A set of strings containg the mark-tipe desired for
%       the line corresponding to each line from model_names.
%   plot_options.u_sel: A set of strings of the same dimension as model_names
%       indicating the name (as declared in dynare) of the exogenous variable whose shock 
%       generates the responses to be ploted.
%   plot_options.v_sel: A set of strings containing the name of the
%       variables (as declared in dynare) whose responses will be displayed.
%       The variables toe be plated will be the same for all the alements in
%       model_name. 
%  plot_options.v_adj: A line vector with the same number of elements as
%       strings in v_sel, indicating the scale adjustment that will be
%       applied for each variable in v_sel.
%  plot_options.v_do_cusum: A line vector with the same number of elements as
%       strings in v_sel, indicating with a 1 if the cummulative impulse
%       response for that variable is required and a 0 if not.
%  plot_options.v_div_ss: A line vector with the same number of elements as
%       strings in v_sel, indicating with a 1 if that variable should be
%       divided by its steady state value and a 0 if not.
%  plot_options.horizon: a numebr indicating the horizon of teh responses to
%       be included in the graphs.
%  plot_options.n_col: a number indicating the number of columns in each the graph.
%  plot_options.n_row: a number indicating the number of rows in each the graph.
%  plot_options.grid: equal to 1 if a grid in each graph is desired, 0 otherwise.
%  plot_options.latex: equal to 1 if a latex format is wanted for titles, 0
%       otherwise. If 1 is selected, the next two have to be included.
%  plot_options.u_name (OPTIONAL): A set of strings of the dimension as u_sel 
%       with the LaTex name for each shock to be included in the title of each graph. 
%       If not included, the name in u_sel is used in the title.
%  plot_options.v_name (OPTIONAL): A set of strings of the dimension as v_sel 
%       with the LaTex name for each variable to be included in the title of each graph. 
%       If not included, the name in v_sel is used in the title.
%
% The function uses the sub-funtion loc defined at the end.

tt=plot_options.horizon;
ir_all=nan(tt,length(plot_options.v_sel));

figure('Name',char(plot_options.u_sel(1)),'units','normalized','outerposition',[0 0 1 1])
set(gcf,'WindowStyle','docked');

for jj=1:length(plot_options.v_sel)
    %subplot(plot_options.n_row,plot_options.n_col,jj);
    subplot_tight(plot_options.n_row,plot_options.n_col,jj,[0.1, 0.06]);
    va = plot_options.v_sel{jj};
    for ii=1:length(plot_options.model_names)

        if isfile([plot_options.model_names{ii} '.mat'])
            load([plot_options.model_names{ii} '.mat'])
            sh = plot_options.u_sel{ii};
            if isempty(loc(M_.endo_names,va))==0
                if plot_options.v_div_ss(jj)==1
                    div_ss=abs(oo_.steady_state(loc(M_.endo_names,va)));
                elseif plot_options.v_div_ss(jj)==2
                    div_ss=abs(exp(oo_.steady_state(loc(M_.endo_names,'yH'))));
                else
                    div_ss=1;
                end
                if isfield(plot_options,'unit_irf')==1
                    if plot_options.unit_irf>0
                        irf_adj=plot_options.unit_irf/(M_.Sigma_e(loc(M_.exo_names,['u_' sh]),loc(M_.exo_names,['u_' sh]))^.5);
                    else
                        irf_adj=1;
                    end
                else
                    irf_adj=1;
                end

                if plot_options.v_do_cusum(jj)==0
                    ir = getfield(oo_.irfs,[va '_u_' sh])...
                        *plot_options.v_adj(jj)/div_ss*irf_adj;
                elseif plot_options.v_do_cusum(jj)==1
                    ir = cumsum(getfield(oo_.irfs,[va '_u_' sh]))...
                        *plot_options.v_adj(jj)/div_ss*irf_adj;
                end

                if max(abs(ir))<1e-12
                    ir=ir*0;
                end
                ir_all(:,ii)=ir(1:tt);
            else
                ir_all(:,ii)=0;
            end
        else
            ir_all(:,ii)=nan;
        end
    end
    for ii=1:length(plot_options.model_names)
        hold on;
        plot(1:tt,ir_all(:,ii),[plot_options.marks{ii}],'LineWidth',2);
        
    end
    if plot_options.grid==1
        grid on;
    end
    if isfield(plot_options,'include_zero')==1
        if plot_options.include_zero==1
            plot(1:tt,0*ir_all(:,ii),'.c','LineWidth',0.5);
        end
    end
    if plot_options.latex==1
        va_name=plot_options.v_name{jj};
        sh_name=plot_options.u_name{1};
    else
        va_name=va;
        sh_name=sh;
    end
    title(['$' sh_name ' \Rightarrow ' va_name '$'],'interpreter','latex','FontSize',14);
    set(gca,'FontSize',12,'FontName','Times')
    xlim([1 tt]);
    hold off;
    %     if jj==length(plot_options.v_sel)
    %         expr=['legend('];
    %         for ii=1:length(plot_options.model_names)
    %             if ii==length(plot_options.model_names)
    %                 expr=[expr plot_options.legend_name{ii} ];
    %             else
    %                 expr=[expr plot_options.legend_name{ii} ','];
    %             end
    %         end
    %         expr=[expr ');'];
    %         eval(expr)
    %     end
end


% MATLAB function: loc.m                        July 22, 1992
%        loc(mstring,'sstring') returns the number of the row of
%        mstring that has the same non-blanck characters as
%        sstring. mstring is a matrix of characters. Each of its
%        rows corresponds to a "name". sstring is a character
%        string. It is the "name" we are looking for in mstring.
%        Note that sstring must be placed in between single
%        quotation marks.

function [x] = loc(mstring,sstring,switchmod)

if ischar(mstring)==0
    mstring=char(mstring);
end

[rm,cm]=size(mstring);
cs=max(size(sstring));

% If necessary, add blanck columns to sstring so it will have the
%  same number of columns as mstring.
if cm>cs;
    nblancks=cm-cs;
    for i=1:nblancks
        sstring=[sstring,' '];
    end
end

if(cm~=max(size(sstring)))
    disp(['problem with padding ',sstring])
    disp('The character string might be longer than name list')
    mstring
    %return
    pause
end

x=[];
for r=1:rm;
    if(length(find(mstring(r,:)==sstring))==cm)
%     if(strcmp(mstring(r,:),sstring)==cm)
        x=r;
    end
end

if(x==0)
    if(~exist('switchmod')); disp(['Could not find ',sstring]); end
end

function vargout=subplot_tight(m, n, p, margins, varargin)
%% subplot_tight
% A subplot function substitude with margins user tunabble parameter.
%
%% Syntax
%  h=subplot_tight(m, n, p);
%  h=subplot_tight(m, n, p, margins);
%  h=subplot_tight(m, n, p, margins, subplotArgs...);
%
%% Description
% Our goal is to grant the user the ability to define the margins between neighbouring
%  subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the
%  margins between subplots can reach 40% of figure area, which is pretty lavish. While at
%  the begining the function was implememnted as wrapper function for Matlab function
%  subplot, it was modified due to axes del;etion resulting from what Matlab subplot
%  detected as overlapping. Therefore, the current implmenetation makes no use of Matlab
%  subplot function, using axes instead. This can be problematic, as axis and subplot
%  parameters are quie different. Set isWrapper to "True" to return to wrapper mode, which
%  fully supports subplot format.
%
%% Input arguments (defaults exist):
%   margins- two elements vector [vertical,horizontal] defining the margins between
%        neighbouring axes. Default value is 0.04
%
%% Output arguments
%   same as subplot- none, or axes handle according to function call.
%
%% Issues & Comments
%  - Note that if additional elements are used in order to be passed to subplot, margins
%     parameter must be defined. For default margins value use empty element- [].
%  - 
%
%% Example
% close all;
% img=imread('peppers.png');
% figSubplotH=figure('Name', 'subplot');
% figSubplotTightH=figure('Name', 'subplot_tight');
% nElems=17;
% subplotRows=ceil(sqrt(nElems)-1);
% subplotRows=max(1, subplotRows);
% subplotCols=ceil(nElems/subplotRows);
% for iElem=1:nElems
%    figure(figSubplotH);
%    subplot(subplotRows, subplotCols, iElem);
%    imshow(img);
%    figure(figSubplotTightH);
%    subplot_tight(subplotRows, subplotCols, iElem, [0.0001]);
%    imshow(img);
% end
%
%% See also
%  - subplot
%
%% Revision history
% First version: Nikolay S. 2011-03-29.
% Last update:   Nikolay S. 2012-05-24.
%
% *List of Changes:*
% 2012-05-24
%  Non wrapping mode (based on axes command) added, to deal with an issue of disappearing
%     subplots occuring with massive axes.
%% Default params
isWrapper=false;
if (nargin<4) || isempty(margins)
    margins=[0.04,0.04]; % default margins value- 4% of figure
end
if length(margins)==1
    margins(2)=margins;
end
%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  
height=(1-(m+1)*margins(1))/m; % single subplot height
width=(1-(n+1)*margins(2))/n;  % single subplot width
% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   
merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width
merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
pos=[merged_left, merged_bottom, merged_width, merged_height];
if isWrapper
   h=subplot(m, n, p, varargin{:}, 'Units', 'Normalized', 'Position', pos);
else
   h=axes('Position', pos, varargin{:});
end
if nargout==1
   vargout=h;
end