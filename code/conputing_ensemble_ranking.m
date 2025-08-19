function  five = conputing_ensemble_ranking( predict_score_matrix3,interaction_matrix,n,A,new_array)
        Sco3= predict_score_matrix3;
        final_score=Sco3( interaction_matrix==0);
        
        for i=1:n
            q=A(new_array(1,i),1);
            w=A(new_array(1,i),2);
            s_score=Sco3(q,w);
            T=[s_score;final_score];                                                %sort the score
            index3=sort(T,'descend');
            index3_1=find(index3==s_score);                                          %find the sort of changed variate
            if length(index3_1)~=1
                index3_2=mean(index3_1);                                              %decide the number of the same value
            else
                index3_2=index3_1;                                                    %get the average rank between same value
            end
            five(i)=index3_2;
        end

end

