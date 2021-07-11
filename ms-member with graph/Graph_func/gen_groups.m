function [Gk, member_out] = gen_groups(member_in, ind, filter)

    member_out=member_in;
    member_out.G=zeros(1,member_out.M);
    Adj = genere_groups_graph(member_out.m(:,ind), member_out.P(:,:,ind), filter.epsilonx, filter.epsilonv);
    Gk= connected_components(Adj);
    for i=1:length(Gk)
        Gk(i).index=ind(Gk(i).index);
        member_out.G(:,Gk(i).index)=i;
        Gk(i).center=mean(member_in.m(:,Gk(i).index),2);
        Gk(i).P=mean(member_in.P(:,:,Gk(i).index),3);
    end


end