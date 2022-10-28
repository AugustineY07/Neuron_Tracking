reference_4day.mat contains the reference matching pairs by include clusters that passes the KW significance test, has simScore > 1, and has no duplicated pair.  
Each cell col is a shank from 1-4  
Each cell row is the matching between two consecutive days(2019.11.21, 2019.11.22, 2019.12.03, 2019.12.13)  
  
Within each cell:  
column 1: idx from the second day in clusters that passes KW test  
column 2: idx from the first day in clusters that passes KW test  
column 3: simScore(ranked from high to low)  
column 4: rank of the simScore in col 3 among all its simScores  
column 5: the cluster number in the second day  
column 6: the cluster number in the first day  
