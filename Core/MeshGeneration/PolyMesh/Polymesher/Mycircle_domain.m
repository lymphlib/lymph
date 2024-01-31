%%

% [~,~,~,~,P1] = PolyMesher(@CircleInclusion,1000,100);
% [~,~,~,~,P2] = PolyMesher(@Circle,100,100);
% p = [P1;P2];

[Node,Element,Supp,Load,P] = PolyMesher(@Rectangle,10,100);

%%
h  = 0.2;
rho = 2;
p = [];
for teta = 0 : h : 2*pi-h
        x = rho*cos(teta) ;
        y = rho*sin(teta) ;
        p = [p; x y];
 end

hold on
 plot(p(:,1),p(:,2),'b-');



%
clear; clc;
h  = 0.025;
rho = [0.5 + h/2 : h : 1-h/2];
p = [];
for i = 1 : length(rho)
    for teta = 0 : h : 2*pi-h
        x = rho(i)*cos(teta) ;
        y = rho(i)*sin(teta) ;
        p = [p; x y];
    end
end


% plot(p(:,1),p(:,2),'*');
[Node,Element,Supp,Load,P] = PolyMesher(@CircleDomain,100,100, p);

Mesh3.Node = Node;
Mesh3.Element = Element;
Mesh3.Supp = Supp;
Mesh3.Load = Load;
Mesh3.P = P;
save('Mesh3_quad.mat');





