function[l_p_S,h_p_S,med_temp]=Show_LP_HP_SS(l_p,h_p,SS_T,CIR,CIR_S)
h_pH=h_p';
h_pH(3,:)=1;
l_pH=l_p';
l_pH(3,:)=1;

[H_IR2S]=homography_ndlt(CIR,CIR_S);

% Transform the corners' points & lesion boundary location to IR using
% "H_W2IR"

l_p_S= H_IR2S * l_pH;
h_p_S= H_IR2S * h_pH;

for i=1:length(l_p_S)
l_p_S(:,i)=l_p_S(:,i)/l_p_S(3,i);
h_p_S(:,i)=h_p_S(:,i)/h_p_S(3,i);
end


h_p_S=h_p_S';
l_p_S=l_p_S';

% compute the variation image by subtracting the median value 
med_temp= median(median(SS_T));
Var_SS= SS_T- med_temp;


figure;
imagesc(Var_SS,[-1 max(max(Var_SS))]);
hold on;
scatter(l_p_S(:,1), l_p_S(:,2), 10, 'b');
hold on;
scatter(h_p_S(:,1), h_p_S(:,2), 10, 'm','MarkerFaceColor',[1,1,1]);
colormap('jet');
colorbar;

