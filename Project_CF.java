import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Arrays;



public class Project_CF {
	static String UNICODE="UTF-8";
	static String Dataset = "13YM10";
	static int scaleMax = 13;
	static String[] SubData = {"Story", "Acting", "Visual", "Directory"};
	public static void main(String[] args) throws Exception {
		
		for(int e = 0; e < 4; e++) 
		{
			int[][] useritemmatris=ReadIntArrayFromFile(Dataset + SubData[e] + ".txt"); // 13YM20Directory.txt , 13YM20Story.txt , 13YM20Visual.txt  ,  13YM20Acting.txt
			
			double[] mean=Mean(useritemmatris);
			double[][] deviationmean=DeviationFromMean(useritemmatris,mean);
			double[] std=Std(deviationmean);
			double[][] zscore=Zscore(deviationmean,std);
			
			int[] neighbour={200};
			int[] common_rate_limit={0};
			
			String[] similarity_type = {"Euclidean", "Manhattan", "Chebyshev", "Pearson"};
			for(int k = 0; k < 4; k++)
			{
			//String similarity_type="Chebyshev";
				for(int i=0;i<neighbour.length;i++){
					for(int j=0;j<common_rate_limit.length;j++){
					String error=FindErrorOfSystem(useritemmatris,zscore,mean,std,neighbour[i],common_rate_limit[j],similarity_type[k], scaleMax, SubData[e]);
					String mae=error.split(",")[0];
					String rmse=error.split(",")[1];
					System.out.println("Data: " + Dataset + SubData[e] + ", Similarity Type = " + similarity_type[k] + ", Neighbour:"+neighbour[i]+"\nMAE: "+mae.substring(0, 6)+" RMSE: "+rmse.substring(0, 6)+"\n");
					}
				}
			}
			//ArrayList<Integer> topnitemindices=PredictionWithZscoreTopN(5,1,neighboursindeces,zscore,weights,mean,std);
			//WritelnforArray(topnitemindices);
		}
		System.out.println("Islem basariyla tamamlandi!");
	}
	
	private static String FindErrorOfSystem(int[][] useritemmatris,double[][] zscore, double[] mean, double[] std, int neighbourcount,int common_rate_limit,String similarity_type, int scaleMax, String SubData) throws IOException {
		
		String file_name = similarity_type + "_" + Dataset + "_" + SubData + ".txt";
		
		ArrayList<String> Prediction_result= new ArrayList<String>();
		String error="";
		int user_count=useritemmatris.length;
		int item_count=useritemmatris[0].length;
	
		double mae=0;
		int mae_count=0;
		double rmse=0;
		int rmse_count=0;
		
		for(int i=0;i<user_count;i++){
			int active_user=i;
			String prediction_sonuc="";
			HashMap<Integer,Double> weights=ComputeSimilarityforActiveUser(useritemmatris,zscore,active_user,common_rate_limit,similarity_type);
			ArrayList<Integer> neighboursindeces=FindNeighboursforActiveUser(weights,neighbourcount,active_user);
			for(int j=0;j<item_count;j++){
				double prediction;
				int active_item=j;	
				if(useritemmatris[active_user][active_item]!=99.0){
					double original_value=useritemmatris[active_user][active_item];
					prediction=PredictionWithZscoreForAnItem(active_user,active_item,neighboursindeces,zscore,weights,mean,std);
					if(prediction>scaleMax)
						prediction=scaleMax;
					if(prediction<1)
						prediction=1;
					mae+=Math.abs(original_value-prediction);
					mae_count++;
					rmse+=Math.pow(original_value-prediction, 2);
					rmse_count++;

				}
				else prediction=99;
				
				prediction_sonuc+=prediction;
				
				if(!(j==item_count-1))
				prediction_sonuc+=",";
				
			}
			Prediction_result.add(prediction_sonuc);
		}
		
		WriteListToFile(Prediction_result,file_name,false);
		
		mae=mae/mae_count;
		rmse=Math.sqrt(rmse/rmse_count);
		error=mae+","+rmse;

		return error;
	}




	private static ArrayList<Integer> PredictionWithZscoreTopN(int activeuser,int topN,ArrayList<Integer> nindeces, double[][] zscore,HashMap<Integer,Double> weights, double[] mean, double[] std) {
		ArrayList<Integer> emptyitemindeces = new ArrayList<Integer>();
		ArrayList<Integer> recommendindeces = new ArrayList<Integer>();
		HashMap<Integer, Double> predictionvalues=new HashMap<Integer, Double>();
		int itemcount=zscore[0].length;
		for(int i=0;i<itemcount;i++){
			if((zscore[activeuser][i]==99.0))
				emptyitemindeces.add(i);
		}
		
		for(int i=0;i<emptyitemindeces.size();i++){
			int emptyitemindex=emptyitemindeces.get(i);
			double prediction=PredictionWithZscoreForAnItem(activeuser,emptyitemindex,nindeces,zscore,weights,mean,std);
			predictionvalues.put(emptyitemindex, prediction);
		}
		
		List<Entry<Integer, Double>> sorteditems=entriesSortedByValues(predictionvalues);
		
		for(int i=0;i<topN;i++){
			recommendindeces.add(i, sorteditems.get(i).getKey());
		}
		
		return recommendindeces;
	}


	public static double  PredictionWithZscoreForAnItem(int activeuser,int activeitem,ArrayList<Integer> nindeces,double[][] zscore,HashMap<Integer,Double> weights,double[] mean,double[] std ){
		double sum=0;
		double sumWeight=0;
		double prediction=0;
		
		if(nindeces.size()==0)
			prediction=mean[activeuser];	
		else{
		for(int i=0;i<nindeces.size();i++){
			int nindex=nindeces.get(i);
			double wau=weights.get(nindex);
			double zuq=zscore[nindex][activeitem];
			if(!(zuq==99.0)){
				sumWeight+=wau;
				sum+=wau*zuq;
			}
		}
		
		if(sum!=0)
		prediction=mean[activeuser]+std[activeuser]*sum/sumWeight;
		else
		prediction=mean[activeuser];
		
		}
		
		return prediction;
    }
	
	public static double  PredictionWithDeviaionForAnItem(int activeuser,int activeitem,ArrayList<Integer> nindeces,double[][] devmean,double[] weights,double[] mean,double[] std ){
		double sum=0;
		double sumWeight=0;
		double prediction=0;
		for(int i=0;i<nindeces.size();i++){
			int nindex=nindeces.get(i);
			double wau=weights[nindex];
			double dev=devmean[nindex][activeitem];
			if(!(dev==99.0)){
				sumWeight+=wau;
				sum+=wau*dev;
			}
		}
		prediction=mean[activeuser]+sum/sumWeight;
		
		return prediction;
    }
	
	
	public static double  PredictionWithNoNormalForAnItem(int activeuser,int activeitem,ArrayList<Integer> nindeces,int[][] useritemmatris,double[] weights){
		double sum=0;
		double sumWeight=0;
		double prediction=0;
		for(int i=0;i<nindeces.size();i++){
			int nindex=nindeces.get(i);
			double wau=weights[nindex];
			double rating=useritemmatris[nindex][activeitem];
			if(!(rating==99.0)){
				sumWeight+=wau;
				sum+=wau*rating;
			}
		}
		prediction=sum/sumWeight;
		
		return prediction;
    }
	
	public static ArrayList<Integer> FindNeighboursforActiveUser(HashMap<Integer,Double> weights,int neighbourcount,int activeuser ){
		ArrayList<Integer> neighbourindeces = new ArrayList<Integer>();
		
		if(weights.size()!=0){
			List<Entry<Integer, Double>> sortedweights=entriesSortedByValues(weights);
			
			if(weights.size()<neighbourcount)
				neighbourcount=weights.size();
			
			for(int i=0;i<neighbourcount;i++)
				neighbourindeces.add(i, sortedweights.get(i).getKey());
		}
		
		return neighbourindeces;
    }
	
	public static HashMap<Integer,Double> ComputeSimilarityforActiveUser(int[][] useritemmatris,double[][] zscore,int activeuserindex,int common_rate_limit,String similarity_type){
		HashMap<Integer,Double> weights=new HashMap<Integer,Double>();
		int usercount=zscore.length;
	
		for(int j=0;j<usercount;j++){
			if(j!=activeuserindex)
				if(ComputeSimilarityforTwoUser(useritemmatris,zscore,activeuserindex,j,common_rate_limit,similarity_type)!=-1000)
				weights.put(j,ComputeSimilarityforTwoUser(useritemmatris,zscore,activeuserindex,j,common_rate_limit,similarity_type));
		}
		
		return weights;
    }
	
	public static double ComputeSimilarityforTwoUser(int[][] useritemmatris,double[][] zscore,int firstuserindex,int seconduserindex,int common_rate_limit,String similarity_type){
		double weight=0.0;

		if(similarity_type.equals("Pearson"))
			weight=ComputeSimilarityWithPearson(zscore,firstuserindex,seconduserindex,common_rate_limit);
		if(similarity_type.equals("Euclidean")||similarity_type.equals("Manhattan")||similarity_type.equals("Chebyshev"))
			weight=ComputeSimilaritywithDistanceMeasure(zscore,useritemmatris,firstuserindex,seconduserindex,similarity_type);
		
		return weight;
    }

	public static double ComputeSimilarityWithPearson(double[][] zscore,int activeuserindex,int otheruserindex,int common_rate_limit){
		int itemcount=zscore[0].length;
		double weight=0;
		int count=0;

		for(int j=0;j<itemcount;j++){
			if(!(zscore[activeuserindex][j]==99.0) & !(zscore[otheruserindex][j]==99.0)){
				weight+=zscore[activeuserindex][j]*zscore[otheruserindex][j];
				count++;
			}
		}

		if(count<common_rate_limit)
			weight= -1000;

		return weight;
	}




	public static double ComputeSimilaritywithDistanceMeasure(double[][] zscore,int[][] useritemmatris,int firstuserindex,int seconduserindex,String distance_measure){
		double distance=0.0;

		if(distance_measure.equals("Euclidean"))
			distance=ComputeDistancewithEuclidean(zscore,useritemmatris,firstuserindex,seconduserindex);
		else if(distance_measure.equals("Manhattan"))
			distance=ComputeDistancewithManhattan(zscore,useritemmatris,firstuserindex,seconduserindex);
		else if(distance_measure.equals("Chebyshev"))
			distance=ComputeDistancewithChebyshev(zscore,useritemmatris,firstuserindex,seconduserindex);

		double similarity_weight=1/(1+distance);
		return similarity_weight;
	}



	public static double ComputeDistancewithEuclidean(double[][] zscore,int[][] useritemmatris,int firstuserindex,int seconduserindex){
		int itemcount=zscore[0].length;
		int total=0;
		for(int j=0;j<itemcount;j++){
			if(!(zscore[firstuserindex][j]==99.0) & !(zscore[seconduserindex][j]==99.0))
				total+=Math.pow(Math.abs(zscore[firstuserindex][j]-zscore[seconduserindex][j]),2);
		}
		double distance=Math.sqrt(total);

		return distance;
	}

	public static double ComputeDistancewithManhattan(double[][] zscore,int[][] useritemmatris,int firstuserindex,int seconduserindex){
		int itemcount=zscore[0].length;
		double total=0.0;
		for(int j=0;j<itemcount;j++){
			if(!(zscore[firstuserindex][j]==99.0) & !(zscore[seconduserindex][j]==99.0))
				total+=Math.abs(zscore[firstuserindex][j]-zscore[seconduserindex][j]);
		}

		return total;
	}

	public static double ComputeDistancewithChebyshev(double[][] zscore,int[][] useritemmatris,int firstuserindex,int seconduserindex){
		int itemcount=zscore[0].length;
		double max=-1;
		for(int j=0;j<itemcount;j++){
			if(!(zscore[firstuserindex][j]==99.0) & !(zscore[seconduserindex][j]==99.0)) {
				double value=Math.abs(zscore[firstuserindex][j] - zscore[seconduserindex][j]);
				if(value>max)
					max=value;
			}
		}
		return max;
	}


	public static double[][] Zscore(double[][] deviationmean,double [] std){
		int usercount=deviationmean.length;
		int itemcount=deviationmean[0].length;
		double[][] zscore=new double[usercount][itemcount];
		for(int i=0;i<usercount;i++){
			for(int j=0;j<itemcount;j++){
				if(!(deviationmean[i][j]==99.0)){
					zscore[i][j]=deviationmean[i][j]/std[i];
				}
				else{
					zscore[i][j]=99.0;
				}
			}
		}
		return zscore;
    }
	
	public static double[] Mean(int[][] useritemmatris){
			int usercount=useritemmatris.length;
			int itemcount=useritemmatris[0].length;
			double[] mean=new double[usercount];
		
		for(int i=0;i<usercount;i++){
			int toplam=0;
			int count=0;
			for(int j=0;j<itemcount;j++){
				if(!(useritemmatris[i][j]==99)){
					toplam+=useritemmatris[i][j];
					count++;
				}
			}
			mean[i]=(double)toplam/count;
		}
		return mean;
    }
	
	public static double[] Std(double[][] deviationfrommean){
		int usercount=deviationfrommean.length;
		int itemcount=deviationfrommean[0].length;
		double[] std=new double[usercount];
		for(int i=0;i<usercount;i++){
			double toplam=0;
			int count=0;
			for(int j=0;j<itemcount;j++){
				if(!(deviationfrommean[i][j]==99.0)){
					toplam+=deviationfrommean[i][j]*deviationfrommean[i][j];
					count++;
				}
			}
			std[i]=Math.sqrt((double)toplam/(count-1));
		}
		return std;
    }
	
	public static double[][] DeviationFromMean(int[][] useritemmatris,double [] mean){
		int usercount=useritemmatris.length;
		int itemcount=useritemmatris[0].length;
		double[][] devmean=new double[usercount][itemcount];
		for(int i=0;i<usercount;i++){
			for(int j=0;j<itemcount;j++){
				if(!(useritemmatris[i][j]==99)){
					devmean[i][j]=useritemmatris[i][j]-mean[i];
				}
				else{
					devmean[i][j]=99;
				}
			}
		}
		return devmean;
    }
	
	public static void WritelnforArray(double[] data){
	    for(int i=0;i<data.length;i++)
	    	System.out.println(data[i]);  
	}
	public static void WritelnforArray(int[] data){
	    for(int i=0;i<data.length;i++)
	    	System.out.println(data[i]);  
	}
	public static void WritelnforArray(ArrayList<Integer> data){
	    for(int i=0;i<data.size();i++)
	    	System.out.println(data.get(i));  
	}
	public static void WriteDoubleMultiArray(double[][] veri){
		for(int i=0;i<veri.length;i++){
			for(int j=0;j<veri[0].length;j++){
		    	System.out.print(veri[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	public static int[][] ReadIntArrayFromFile(String input_file) throws IOException {	
		ArrayList<String> datalist=ReadFileAsList(input_file);
		String data = datalist.get(0);
		String[] data_string=data.split(",");
		int[][] datarows=new int[datalist.size()][data_string.length] ;
		for(int i=0;i<datalist.size();i++){
			data = datalist.get(i);
			data_string=data.split(",");
			for(int k=0;k<data_string.length;k++){
				datarows[i][k]=(int) atof(data_string[k]);
			}
		}
		return datarows;
	}
	
	public static ArrayList<String> ReadFileAsList(String filename) throws IOException //ok
    {
		ArrayList<String> list = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filename),UNICODE));
		String line;

        while ((line = br.readLine()) != null) 
        	list.add(line);
        
        br.close();
        return list;
    }
	public static double atof(String s){
		return Double.valueOf(s).doubleValue();
	}
	public static <K,V extends Comparable<? super V>> List<Entry<K, V>> entriesSortedByValues(Map<K,V> map) {

		List<Entry<K,V>> sortedEntries = new ArrayList<Entry<K,V>>(map.entrySet());

		Collections.sort(sortedEntries,new Comparator<Entry<K,V>>() {
			public int compare(Entry<K,V> e1, Entry<K,V> e2) {
            return e2.getValue().compareTo(e1.getValue());
        }
		});
		
		return sortedEntries;
	}
	
	public static void WriteListToFile(ArrayList<String> list,String dosya_yolu,boolean ustuneyaz) throws IOException //ok
    {  
        try {
            FileOutputStream fos = new FileOutputStream(dosya_yolu,ustuneyaz);
            Writer out = new OutputStreamWriter(fos, UNICODE);
            
            for(int i=0;i<list.size();i++){
    			out.append(list.get(i));
    			
    			if(!(i==list.size()-1))
    				out.append("\n");
    		}
            
            out.close();
        } 
        catch (IOException e) {
            e.printStackTrace();
        }
        return;  
    }
	
	public static void DivideDataForCrossValidation(String datapath, int crossnumber,String outputpath) throws IOException {
		ArrayList<String> folders=ListFolderNamesinDirectory(datapath);
		for(int i=1;i<=crossnumber;i++){ // create empty folders
			CreateNewFolder(outputpath+"/"+i);
			CreateNewFolder(outputpath+"/"+i+"/"+"training");
			CreateNewFolder(outputpath+"/"+i+"/"+"test");
			for(int j=0;j<folders.size();j++){
				CreateNewFolder(outputpath+"/"+i+"/"+"training"+"/"+folders.get(j));
				CreateNewFolder(outputpath+"/"+i+"/"+"test"+"/"+"/"+folders.get(j));
			}
		}
		
		for(int k=0;k<folders.size();k++){
			ArrayList<String> files=ListFileNamesinDirectory(datapath+"/"+folders.get(k));
			Collections.shuffle(files);
			
			int filesize =files.size();
			int partsize=filesize/crossnumber;
			ArrayList<ArrayList<String>> splittedfiles = new ArrayList<ArrayList<String>>();
			
			for(int i=0;i<crossnumber;i++){
				List<String> subfiles=files.subList(i*partsize, (i+1)*partsize);
				ArrayList<String> templist=new ArrayList<String>();
				templist.addAll(subfiles);
				splittedfiles.add(i,templist);
			}
			ArrayList<ArrayList<String>> testsets=splittedfiles;
			
			ArrayList<ArrayList<String>> trainingsets=new ArrayList<ArrayList<String>>();
			
			for(int i=0;i<splittedfiles.size();i++){
				ArrayList<String> trainingset=new ArrayList<String>();
				for(int j=0;j<splittedfiles.size();j++){
					if(i!=j)
					trainingset.addAll(splittedfiles.get(j));
				}
				trainingsets.add(i, trainingset);
			}
			
			for(int i=0;i<testsets.size();i++){
				ArrayList<String> testset=testsets.get(i);
				for(int j=0;j<testset.size();j++){
					CopyFile(datapath+"/"+folders.get(k)+"/"+testset.get(j),outputpath+"/"+(i+1)+"/"+"test"+"/"+folders.get(k)+"/"+testset.get(j));
				}
			}
			
			for(int i=0;i<trainingsets.size();i++){
				ArrayList<String> trainset=trainingsets.get(i);
				for(int j=0;j<trainset.size();j++){
					CopyFile(datapath+"/"+folders.get(k)+"/"+trainset.get(j),outputpath+"/"+(i+1)+"/"+"training"+"/"+folders.get(k)+"/"+trainset.get(j));
				}
			}
			
		}
		
}
	
	public static ArrayList<String> ListFolderNamesinDirectory(String dosyaismi) {
		File folder = new File(dosyaismi);
		File[] listOfFiles = folder.listFiles();
		ArrayList<String> Folders=new ArrayList<String>();

		    for (int i = 0; i < listOfFiles.length; i++) {
		      if (listOfFiles[i].isDirectory())
		    	  Folders.add(listOfFiles[i].getName());
		    }
		return Folders;
	}
	
	public static void CopyFile(String sourcepath, String destinationpath) throws IOException {
		WriteListToFile(ReadFileAsList(sourcepath), destinationpath, false);
	}
	
	public static void CreateNewFolder(String dosyayolu){
		boolean success = (new File(dosyayolu)).mkdirs();
        if (!success)  System.out.println("Dosya Olusturulamadi");
	}
	
	public static ArrayList<String> ListFileNamesinDirectory(String dosyaismi) {
		File folder = new File(dosyaismi);
		File[] listOfFiles = folder.listFiles();
		ArrayList<String> Files=new ArrayList<String>();

		    for (int i = 0; i < listOfFiles.length; i++) {
		      if (listOfFiles[i].isFile())
		    	  Files.add(listOfFiles[i].getName());
		    }
		return Files;
	}
	
	
	
}
