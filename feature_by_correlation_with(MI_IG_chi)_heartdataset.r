#data<-read.csv("/home/pratik/Downloads/Data/breastcancer.csv",header = T)
data<-read.csv("/home/pratik/Downloads/Data/heart.csv",header = T)
require(infotheo);require(igraph)
require(entropy);require(modMax);require(FSelector)

#-------------------------------------------find mutual information--------------------------------------------#
data_class<-data[,dim(data)[2]]
data_class[which(data_class == 1)]<- "B"
data_class[which(data_class == 0)]<- "M"
response1<-as.numeric(data[,dim(data)[2]])
mut_info1<-c()
for(i in 1:(dim(data)[2]-1)){
  dis_cret<-discretize2d(data[,i],response1,10,2)
  mut<-mi.empirical(dis_cret)
  mut_info1<-c(mut_info1,mut)
}

#-------------------------------------------find chi-square value----------------------------------------------#
chi_value<-chi.squared(target~.,data)

#-------------------------------------------find information gain----------------------------------------------#
information_value<-information.gain(target~.,data)

average_value<-sapply(1:(dim(data)[2]-1),function(i)mean(mut_info1[i],chi_value[i,1],information_value[i,1]))

#---------------------------------------------thresold value---------------------------------------------------#
l<-length(which(average_value > 0.06))
data<-data[,c(which(average_value > 0.06))]

#--------------------------------------spearman correlation coefficient----------------------------------------#
details<-data.frame()
for(j in 1:(l-1)){
  m<-seq((j+1),l,1)
  for(k in m){
    cor_value<-cor(data[,j],data[,k],method = "spearman")
    value<-data.frame(feature1 = j,feature2 = k,spearman_cor = cor_value)
    details<-rbind(value,details)
    #   spearm_cor<-c(spearm_cor,cor_value)
  }
}

#--------------------------------------------ordered vertex---------------------------------------------------#
sorted_feature<-data.matrix(details[order(details$spearman_cor,decreasing = T),])
n<-dim(sorted_feature)[1]

feature_vertex<-c()
for(t in 1:n){
  feature_ver<-as.vector(c(sorted_feature[t,][1],sorted_feature[t,][2]))
  feature_vertex<-c(feature_vertex,feature_ver)
}

#---------------------------------------------graph generation------------------------------------------------#
i<-1
repeat{
  g<-graph(feature_vertex[1:(2*i)],directed = F)
  c<-components(g)
  number_vertex<-c$csize
  number_component<-c$no
  if( number_vertex ==  l && number_component == 1){
    break
  }
  else{
    i<-i+1 
  }
}

plot(g,vertex.color = "lightgreen",edge.color = "blue")

#-----------------------------------------------find cluster--------------------------------------------------#
cluster_1<-cluster_fast_greedy(g, merges = T, modularity = T)
plot(cluster_1,g,vertex.size = 10,vertex.label.cex = 1.5)
cluster_num_1<-length(cluster_1)
comm_1<-communities(cluster_1)
number_cluster_1<-dim(comm_1)

#-------------------------------------------final feature selection-------------------------------------------#
feature_selected<-c()

for(q in 1:number_cluster_1){
  feature_vertices<-comm_1[[q]]
  nf<-length(feature_vertices)
  
  while(nf > 0 ){
    feature_mutual_info<-mut_info1[feature_vertices]
    feature_selection<-feature_vertices[which.max(feature_mutual_info)]
    feature_selected<-c(feature_selection,feature_selected)
    remove_feature<-c(feature_selection,neighbors(g,feature_selection))
    feature_vertices<-setdiff(feature_vertices,remove_feature)
    nf<-length(feature_vertices)
  }
}


data<-cbind(data[,feature_selected],data_class)
suffle_data<-data[sample.int(nrow(data)),]
n_r<-floor(0.8*nrow(data))
training<-suffle_data[1:n_r,]
test<-suffle_data[-(1:n_r),]

#------------------------------------apply algorithm for classification----------------------------------#
require(randomForest)

kfold<-5
r<-floor(nrow(suffle_data)/kfold)
accuracy<-c()

for( i in 1:kfold){
  test_data<- suffle_data[((r*(i-1)+1):(r*i)),]
  training_data<-suffle_data[-((r*(i-1)+1):(r*i)),]
  
  rf<-randomForest(data_class~.,training_data ,importance = TRUE,proximity = TRUE)
  test_pred<-predict(rf,test_data)
  acc<-mean(test_pred == test_data$data_class)
  accuracy<-c(acc,accuracy)
}
accuracy
max(accuracy)

p<-which.max(accuracy)
final_training_data<-suffle_data[-((r*(p-1)+1):(r*p)),]
rf2<-randomForest(data_class ~.,training_data ,importance = TRUE,proximity = TRUE,nPerm =3)


