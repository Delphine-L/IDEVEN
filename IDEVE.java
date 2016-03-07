import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;


public class IDEVE {

	private static Vector<String> family_list = new Vector<String>();
	private static String member;
	private static String fasta_line;
	private static ArrayList<String> othoList = new ArrayList<String>();
	private static String sep = "[|]+";
	private static String line_orth;
	private static String info1[];
	private static String info2[];
	private static String paire_o;
	private static String paire_o_r;
	private static String gene1;
	private static String gene2;
	private static String[] temp_array1;
	private static String[] temp_arrayb;
	private static Vector<String> familyblocs= new Vector<String>();
	private static String[] temp_array;
	
	private static String bloctemps;
	private static String ksfilename; //Ks files inter species
	private static String inpara1ks; //Ks file within species 1
	private static String inpara2ks; //Ks file within species 2
	private static String inparalogous2_bloc; //Blocks file within species 2
	private static String inparalogous1_bloc; //Blocks file within species 1
	private static String paires_blocs; //Blocks file interspecies ,  all syntenic depth
	private static String ortho_blocs; //Blocks file interspecies ,  selected syntenic depth
	private static String tandem1_file ; //Tandem duplications species 1
	private static String tandem2_file ; //Tandem duplications species 2
	private static String familyfile ; // Family file
	private static String outputfile;//output
	
	
	private static  String analysed_pair = new String();
	private static  String analysed_pair_r = new String();
//	String[] analysed_bloc;
	
	private static  int t=0;
	private static  int c=0;
	private static  int fam_count =0 ;
	private static  int wgd_count =0 ;
	private static  String type = "";
	private static  double KS=0.0;
	private static  String analysed_bloc;
	private static  String neighbour_pairs;
	private static  String[] neighbour_pairs_list;
	private static  int nbr_synt_bloc;
	private static  int nbr_tot_bloc;
	private static  Double neighb_ks;
	private static  Double med_ks =0.0;
	private static  String genea;
	private static  String geneb;
	private static  String[] wgdlist;
	private static  String analysed_wgd_pair;
	private static  String analysed_wgd_pair_r;
	private static  String[] gene_array = new String[2];
	private static  String name1;
	private static  double KS2=0.0;
	private static  String analysed_bloc2;
	private static  String neighbour_pairs2;
	private static  String[] neighbour_pairs_list2;
	private static  int nbr_synt_bloc2;
	private static  double neighb_ks2;
	private static  double med_ks2 =0.0;
	
	
	private static HashMap<String, String> genelist =new HashMap<String, String> (); // contient  nom et info genes inter especes
	private static HashMap<String, String> pairInf = new HashMap<String, String>(); //contient paires et bloc inter especes 
	private static HashMap<String, String> syntlist =new HashMap<String, String> (); //contient un gene, et liste gene synt inter especes 
	private static HashMap<String, String> bloclist =new HashMap<String, String> (); //contient un bloc, et liste genes compris inter especes 
	
	
	private static HashMap<String, String> para1_genelist =new HashMap<String, String>(); // contient  nom et info genes
	private static HashMap<String, String> para1_pairInf = new HashMap<String, String>(); //contient paires et bloc
	private static HashMap<String, String> para1_syntlist =new HashMap<String, String>(); //contient un gene, et liste gene synt
	private static HashMap<String, String> para1_bloclist =new HashMap<String, String>(); //contient un bloc, et liste genes compris
		
	private static HashMap<String, String> para2_genelist =new HashMap<String, String>(); // contient  nom et info genes
	private static HashMap<String, String> para2_pairInf = new HashMap<String, String>(); //contient paires et bloc
	private static HashMap<String, String> para2_syntlist =new HashMap<String, String>(); //contient un gene, et liste gene synt
	private static HashMap<String, String> para2_bloclist =new HashMap<String, String>(); //contient un bloc, et liste genes compris
	
	private static HashMap<String, String> ortho_genelist =new HashMap<String, String>(); // contient  nom et info genes
	private static HashMap<String, String> ortho_pairInf = new HashMap<String, String>(); //contient paires et bloc
	private static HashMap<String, String> ortho_syntlist =new HashMap<String, String>(); //contient un gene, et liste gene synt
	private static HashMap<String, String> ortho_bloclist =new HashMap<String, String>(); //contient un bloc, et liste genes compris
	
	private static HashMap<String, Double> pairList = new HashMap<String, Double>(); //contient paires et Ks
	private static HashMap<String, Double> para1List = new HashMap<String, Double>();
	private static HashMap<String, Double> para2List = new HashMap<String, Double>();
	
	public static void main(String[] args) {
		// command line java -jar IDEVE_test.jar all_ks.txt inpara.txt inpara_sorgho.txt inpara_sorgho_bloc.txt inpara_bloc.txt all_blocs.txt family.txt ortho.txt test_result.txt

		BufferedReader paires_ks = null;
		BufferedReader family = null;
//		BufferedReader orthologous = null;
		BufferedReader inparalogous1 = null;
		BufferedReader inparalogous2 = null;
		BufferedWriter resultFile=null;
		BufferedReader orthologous = null;
		BufferedReader tandem1=null;
		BufferedReader tandem2 = null;
		
	    try{
//	    	String ksfilename= "/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe-tomate_all_ks.txt";
//	    	String inpara1ks="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/tomate_ks.txt";
//	    	String inpara2ks="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe_ks.txt";
//	    	String orthoks="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe-tomate_1-3_Ks.txt";
//	    	String inparalogous2_bloc="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe_blocs.txt";
//	    	String inparalogous1_bloc="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/tomate_blocs.txt";
//	    	String paires_blocs="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe-tomate_all_blocs.txt";
//	    	family = new BufferedReader(new FileReader("/Users/delphine/Documents/test intreegreat/Galaxy204-[Concatenate_datasets_on_data_203_and_data_10].fasta"));
//	    	orthologous  = new BufferedReader(new FileReader("/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe-tomate_1-3_blocs.txt"));
//	    	File file= new File("/Users/delphine/Documents/resultats_CoGe/result_caffe_tomate.txt");
//	    	String ortho_blocs="/Users/delphine/Documents/resultats_CoGe/caffe_tomate/caffe-tomate_1-3_blocs.txt";
//	    	File json_MAIZE= new File("/Users/delphine/Dropbox/workspace2/MAIZE-MGDB5b60-locus_tag.json");
//	    	File json_SORBI= new File("/Users/delphine/Dropbox/workspace2/SORBI-phytozome8_0-locus_tag.json");
//	    	//	    	orthologous  = new BufferedReader(ortho_blocs);
	    	ksfilename= args[0]; //Ks files inter species
	    	inpara1ks=args[1]; //Ks file within species 1
	    	inpara2ks=args[2]; //Ks file within species 2
	    	inparalogous2_bloc=args[3]; //Blocks file within species 2
	    	inparalogous1_bloc=args[4]; //Blocks file within species 1
	    	paires_blocs=args[5]; //Blocks file interspecies ,  all syntenic depth
	       	ortho_blocs=args[6]; //Blocks file interspecies ,  selected syntenic depth
	    	tandem1_file=args[7]; //Tandem duplications species 1
	    	tandem2_file=args[8]; //Tandem duplications species 2
	    	familyfile=args[9]; // Family file
	    	outputfile=args[10]; //output
	    	
	    	
	    	
	    	File file= new File(outputfile); //output
	    	orthologous  = new BufferedReader(new FileReader(ortho_blocs));
	    	resultFile = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
	    	resultFile.write("#Nom gene1\tNom Gene2\tevent\tKs\tMean Ks\tsize bloc\n");
	     	

	    	
			//===========================================================================Recup tandem content ===========================================================================//
	    	
	    	
	    	ArrayList<String> tandem1_list = gettandem(tandem1_file); //liste de paires concat�n�es
			ArrayList<String> tandem2_list = gettandem(tandem2_file);
			
	    	//===========================================================================Recup famille ===========================================================================//
	    	
			family = new BufferedReader(new FileReader(familyfile)); // Family file
	    	while ((fasta_line = family.readLine()) != null) { //parcours chaque ligne du ficher
	    		if(fasta_line.length()!=0 && fasta_line.substring(0,1).matches(">")){ //si en tete de fasta
	    			member=fasta_line.substring(1,fasta_line.length()); //recup�ration du nom
	    			member=member.toUpperCase(); // passage en majuscule
	    			family_list.add(member); // ajout � la liste des membres
	    			//System.out.println(member);
	    		}
			}
	    	
	    	//===========================================================================Recup info fichier Synmap ===========================================================================//

	    	
	    	ArrayList<HashMap<String, String>> finallist=parseblocfile(paires_blocs);
	    	genelist = finallist.get(0); // contient  nom et info genes
	    	pairInf = finallist.get(1); //contient paires et bloc
	    	syntlist = finallist.get(2); //contient un gene, et liste gene synt
	    	bloclist = finallist.get(3); //contient un bloc, et liste genes compris

	    	
	    	//===========================================================================Recup info fichier inpara1 ===========================================================================//
	    	

	    	
	    	ArrayList<HashMap<String, String>> finalpara1=parseblocfile(inparalogous1_bloc);
	    	para1_genelist = finalpara1.get(0); // contient  nom et info genes
	    	para1_pairInf = finalpara1.get(1); //contient paires et bloc
	    	para1_syntlist = finalpara1.get(2); //contient un gene, et liste gene synt
	    	para1_bloclist = finalpara1.get(3); //contient un bloc, et liste genes compris

	    	
	    	//===========================================================================Recup info fichier ortho_blocs ===========================================================================//
	    	
	    	
	    	ArrayList<HashMap<String, String>> finalortho=parseblocfile(ortho_blocs);
	    	ortho_genelist = finalortho.get(0); // contient  nom et info genes
	    	ortho_pairInf = finalortho.get(1);  //contient paires et bloc
	    	ortho_syntlist = finalortho.get(2); //contient un gene, et liste gene synt
	    	ortho_bloclist = finalortho.get(3);  //contient un bloc, et liste genes compris
	    	
	    	
	    	//===========================================================================Recup info fichier inpara2 ===========================================================================//
	    		    	
	    	ArrayList<HashMap<String, String>> finalpara2=parseblocfile(inparalogous2_bloc);
	    	para2_genelist = finalpara2.get(0); // contient  nom et info genes
	    	para2_pairInf = finalpara2.get(1);   //contient paires et bloc
	    	para2_syntlist = finalpara2.get(2);  //contient un gene, et liste gene synt
	    	para2_bloclist = finalpara2.get(3);  //contient un bloc, et liste genes compris

	    	
	    	//--------------------------------------------------------------Paires et Ks -----------------------------------------------------------//
	    	
	    	System.out.println("Find KS");
	    	
	    	pairList=getks(ksfilename); //contient paires et Ks Inter especes 
	    	para1List= getks(inpara1ks); //contient paires et Ks espece 1
	    	para2List= getks(inpara2ks); //contient paires et Ks espece 2
	    	
	    	//===========================================================================Completion famille et recup blocs famille ===========================================================================//
	    	
	    	//Vector<String> temp=new Vector<String>();
//	    	String[] temp_array2;
	    	
	    	for(Integer compt = 0; compt < family_list.size(); compt++) //parcourir chaque element de la liste contenant les membres de la famille
				{
				    String gene = family_list.get(compt);
				    gene=gene.trim();
				    gene=gene.toUpperCase();
				    //System.out.println(gene);
				    //temp_array=gene.split("_");
				    
				    //if (temp_array.length>=2){
				    //	temp_array[temp_array.length-1]=temp_array[temp_array.length-1].split("\\.")[0];
				    //	temp_array[temp_array.length-1]=temp_array[temp_array.length-1].replace("T", "P");
				    //    gene = String.join("_", temp_array);
				    //}
				    if(syntlist.containsKey(gene)){ //Regarder si la liste contenant les genes syntenique � un gene contient le g�ne d'int�ret
					    for( String elem : syntlist.get(gene).split(",")) //parcourir la liste des genes synteniques aux genes
					    {
					    	if(!family_list.contains(elem)) //si le gene syntenique n'est pas present dans la famile, ajouter � la famille
			    			{
					    		elem=elem.toUpperCase();
		    					family_list.add(elem);
			    			}
					    }
					    bloctemps = genelist.get(gene).split(sep)[1]; //recup�rer l'identifiant de blocs syntenique
					    if(!familyblocs.contains(bloctemps))
		    			{
					    	familyblocs.add(bloctemps); //ajout du bloc contenant le membre de la famille � la liste des blocs de la famille
		    			}	    	
				    }else {
				    }
				}
	    	// System.out.println(family_list);
	    	
	    	//===========================================================================creation list events ===========================================================================//
	    	

	    	
	    	while ((line_orth = orthologous.readLine()) != null) { //Parcourir les lignes (paires) du ficher contenant les paires filtr�es par profondeur et class�es par blocs 
		    	paire_o= new String();
		    	info1= new String[9]; 
		    	info2= new String[9];
		    	String[] line_array=line_orth.split("\t");
		    	if (line_array.length >= 9 && !line_array[0].startsWith("#") ) //Recup�rer seulement les lignes contenant les paires (pas les entetes etc)
		    	{ 
		    		info1=line_array[1].split(sep); //recuperer les infos du g�ne 1 et les ranger separement dans une liste
		    		info2=line_array[5].split(sep); //recuperer les infos du g�ne 2 et les ranger separement dans une liste
		    		gene1=info1[3].trim(); 
		    		gene1=gene1.toUpperCase();
		    		gene2=info2[3].trim();
		    		gene2=gene2.toUpperCase();
		    		//temp_array1=gene1.split("_");
		    		//temp_arrayb=gene2.split("_");
				    
				    //if (temp_array1.length>=2){
				    //	temp_array1[temp_array1.length-1]=temp_array1[temp_array1.length-1].split("\\.")[0];
				    //	temp_array1[temp_array1.length-1]=temp_array1[temp_array1.length-1].replace("T", "P");
				    //	gene1 = String.join("_", temp_array1);
				    //}
				    
				   // if (temp_arrayb.length>=2){
				   // 	temp_arrayb[temp_arrayb.length-1]=temp_arrayb[temp_arrayb.length-1].split("\\.")[0];
				   // 	temp_arrayb[temp_arrayb.length-1]=temp_arrayb[temp_arrayb.length-1].replace("T", "P");
				   // 	gene2 = String.join("_", temp_arrayb);
				    //}
		    		
		    		paire_o= gene1+","+gene2; //concatenation des noms de deux genes pour former le nom de la paire, sens 1
		    		paire_o_r= gene2+","+gene1; //sens 2
		    		if(!othoList.contains(paire_o)){  //ajout des paires dans la liste des orthologues
		    			othoList.add(paire_o);
		    		}
		    		if(!othoList.contains(paire_o_r)){
		    			othoList.add(paire_o_r);
		    		}
		    	}
	    	}
	    	

	    	//===========================================================================ID events ===========================================================================//

			String name2;
	    	while(fam_count<family_list.size()-1){ //pour chaque membre de la famille
	    		for(t=fam_count+1; t<family_list.size(); t++){ //Pour chaque autre membre de la famille
	    			//System.out.println(family_list.get(fam_count)+"\t"+family_list.get(t));
	    			if(family_list.get(fam_count).contains(".")){ 
	    				name1=family_list.get(fam_count);//.split("\\.")[0];
	    				name1=name1.toUpperCase();
	    			} else {
	    				name1=family_list.get(fam_count);
	    				name1=name1.toUpperCase();
	    			}
	    			if(family_list.get(t).contains(".")){
	    				//System.out.println(family_list.get(t));
	    				//System.out.println(family_list.get(t).split("\\."));
	    				name2=family_list.get(t);//.split("\\.")[0];
	    				name2=name2.toUpperCase();
	    			} else {
	    				name2=family_list.get(t);
	    				name2=name2.toUpperCase();
	    			}
	    			if( genelist.keySet().contains(name1) && genelist.keySet().contains(name2) ) //teste si les deux genes sont pr�sent dans les listes de g�nes
	    			{
		    			type=new String();
		    			analysed_pair= name1+","+name2;
		    			analysed_pair_r= name2+","+name1;
		    			//System.out.println(analysed_pair_r);
		    			if(othoList.contains(analysed_pair) || othoList.contains(analysed_pair_r)){ //teste si la paire de gènes est dans la liste profondeur syntenique
		    				try {
		    					type="ortholog";
		    				} catch(Exception exp) {
		    					System.out.println("erreur dans type:");
		    					exp.printStackTrace();
		    					
		    				}
		    				try{
		    					KS=pairList.get(analysed_pair); // recuperer le Ks de la paire
		    				}catch(NullPointerException e){
		    					KS=-1; // recuperer le Ks de la paire
		    					System.out.println("Erreur : "+analysed_pair);
						}
		    				analysed_bloc=pairInf.get(analysed_pair);
		    				neighbour_pairs=bloclist.get(analysed_bloc); //recuperer les syntelogs de la paire
		    				neighbour_pairs_list=neighbour_pairs.split(sep);//recuperation de la l iste des paires du bloc
		    				nbr_synt_bloc=neighbour_pairs_list.length/2; //calcul de la taille du bloc
		    				med_ks=0.0;
		    				if (neighbour_pairs_list.length>0){
			    				int err_count=0;
			    				for (String p : neighbour_pairs_list){
			    					try{
			    						med_ks=med_ks+pairList.get(p);
			    					}catch (NullPointerException e){
			    						err_count=err_count+1;
			    					}
			    				}
			    				med_ks=med_ks/(neighbour_pairs_list.length);
		    				}else {
		    					med_ks=0.0;
		    				}
		    				if(type!="")
		    				{
		    					resultFile.write(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
		    				    resultFile.flush();
		    				}
		    				genea=analysed_pair.split(",")[0];
		    				geneb=analysed_pair.split(",")[1];
		    				gene_array[0]=genea;
		    				gene_array[1]=geneb;
		    				for (String gene:gene_array){
		    					try{
		    					wgdlist=ortho_syntlist.get(gene).split(",");
			    				//System.out.println(wgdlist);
		    					wgd_count =0;
			    				if (wgdlist.length>1){
			    					while(wgd_count<wgdlist.length-1){
			    						for(c=wgd_count+1; c<wgdlist.length; c++){
			    							analysed_wgd_pair=wgdlist[wgd_count]+","+wgdlist[c];
			    							analysed_wgd_pair_r=wgdlist[c]+","+wgdlist[wgd_count];
			    							type="WGD";
			    							KS2=0;
			    							if(para1_pairInf.containsKey(analysed_wgd_pair) || para1_pairInf.containsKey(analysed_wgd_pair_r)){
			    								try{
			    									if(para1_pairInf.containsKey(analysed_wgd_pair)){
			    										try{
			    											KS2=para1List.get(analysed_wgd_pair); // recuperer le Ks de la paire
			    					    				}catch(NullPointerException e){
			    					    					KS2=-1; // recuperer le Ks de la paire
			    					    				}
			    										try {
			    			    						    analysed_bloc=para1_pairInf.get(analysed_wgd_pair);
			    			    						}catch(NullPointerException e){
			    					    					analysed_bloc=""; // recuperer le Ks de la paire
			    					    				}
			    			    					}else{
			    			    						try{
			    											KS2=para1List.get(analysed_wgd_pair_r); // recuperer le Ks de la paire
			    					    				}catch(NullPointerException e){
			    					    					KS2=-1; // recuperer le Ks de la paire
			    					    				}
			    			    						try {
			    			    						    analysed_bloc=para1_pairInf.get(analysed_wgd_pair_r);
			    			    						}catch(NullPointerException e){
			    					    					analysed_bloc=""; // recuperer le Ks de la paire
			    					    				}
			    			    						
			    			    					}
			    									neighbour_pairs=para1_bloclist.get(analysed_bloc);
			    				    				neighbour_pairs_list=neighbour_pairs.split(sep);
			    				    				nbr_synt_bloc=neighbour_pairs_list.length/2;
			    				    				med_ks2=0.0;
			    			    					if (neighbour_pairs_list.length>0){
			    					    				int err_count=0;
			    					    				for (String p : neighbour_pairs_list){
			    					    					try{
			    					    						med_ks2=med_ks2+para1List.get(p);
			    					    					}catch (NullPointerException e){
			    					    						err_count=err_count+1;
			    					    					}
			    					    				}
			    					    				med_ks2=med_ks2/(neighbour_pairs_list.length);
			    				    				}else {
			    				    					med_ks2=0.0;
			    				    				}
			    								}catch(NullPointerException e){
			    									KS2=-1;
			    			    					nbr_synt_bloc=0;
			    			    					med_ks2=0.0;
			    								}
			    								System.out.println("RES : "+wgdlist[wgd_count]+"\t"+wgdlist[c]+"\t"+KS2);
			    								if(type!="")
			    								{
			    									resultFile.write(wgdlist[wgd_count]+"\t"+wgdlist[c]+"\t"+type+"\t"+KS2+"\t"+med_ks2+"\t"+nbr_synt_bloc+"\n");
			    									    resultFile.flush();
			    								}
			    								//System.out.println(wgdlist[wgd_count]+"\t"+wgdlist[c]+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
			    							}else if(para2_pairInf.containsKey(analysed_wgd_pair) || para2_pairInf.containsKey(analysed_wgd_pair_r)){
			    								try {
			    			    					if(para2_pairInf.containsKey(analysed_wgd_pair)){
			    			    						try{
			    											KS2=para2List.get(analysed_wgd_pair); // recuperer le Ks de la paire
			    					    				}catch(NullPointerException e){
			    					    					KS2=-1; // recuperer le Ks de la paire
			    					    				}
			    			    						try {
			    			    						    analysed_bloc=para2_pairInf.get(analysed_wgd_pair);
			    			    						}catch(NullPointerException e){
			    					    					analysed_bloc=""; // recuperer le Ks de la paire
			    					    				}
			    			    					}else{
			    			    						try{
			    											KS2=para2List.get(analysed_wgd_pair_r); // recuperer le Ks de la paire
			    					    				}catch(NullPointerException e){
			    					    					KS2=-1; // recuperer le Ks de la paire
			    					    				}
			    			    						try{
			    			    						    analysed_bloc=para2_pairInf.get(analysed_wgd_pair_r);
			    			    						}catch(NullPointerException e){
			    					    					analysed_bloc=""; // recuperer le Ks de la paire
			    					    				}
			    			    					}
			    			    					neighbour_pairs=para2_bloclist.get(analysed_bloc);
			    				    				neighbour_pairs_list=neighbour_pairs.split(sep);
			    				    				nbr_synt_bloc=neighbour_pairs_list.length/2;
			    				    				med_ks2=0.0;
			    				    				if (neighbour_pairs_list.length>0){
			    					    				int err_count=0;
			    					    				for (String p : neighbour_pairs_list){
			    					    					try{
			    					    						med_ks2=med_ks2+para2List.get(p);
			    					    					}catch (NullPointerException e){
			    					    						err_count=err_count+1;
			    					    					}
			    					    				}
			    					    				med_ks2=med_ks2/(neighbour_pairs_list.length);
			    				    				}else {
			    				    					med_ks2=0.0;
			    				    				}
			    								}catch(NullPointerException e){
				    								KS2=-1;
				    			    				nbr_synt_bloc=0;
				    			    				med_ks2=0.0;
				    							}
			    								System.out.println("RES : "+ wgdlist[wgd_count]+"\t"+wgdlist[c]+"\t"+KS2);
			    								resultFile.write(wgdlist[wgd_count]+"\t"+wgdlist[c]+"\t"+type+"\t"+KS2+"\t"+med_ks2+"\t"+nbr_synt_bloc+"\n");
			    								resultFile.flush(); //System.out.println(wgdlist[wgd_count]+"\t"+wgdlist[c]+"\t"+type+"\t"+KS2+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
			    							}
			    						}
			    						wgd_count=wgd_count+1;
			    					}
			    				}
		    				}
			    			catch(NullPointerException e){
			    				System.out.println("Null expression error : "+gene);
			    			
		    				}finally{
		    					
		    				}
		    				}
		    				
		    				
		    				
		    			}else if(para1_pairInf.containsKey(analysed_pair) || para1_pairInf.containsKey(analysed_pair_r)){
		    				type="inParalog";
		    				try {
		    					if(para1_pairInf.containsKey(analysed_pair)){
		    						try{
		    							KS=para1List.get(analysed_pair); // recuperer le Ks de la paire
				    				}catch(NullPointerException e){
				    					KS=-1; // recuperer le Ks de la paire
				    				}
			    					analysed_bloc=para1_pairInf.get(analysed_pair);
		    					}else{
		    						try{
		    							KS=para1List.get(analysed_pair_r); // recuperer le Ks de la paire
				    				}catch(NullPointerException e){
				    					KS=-1; // recuperer le Ks de la paire
				    				}
			    					analysed_bloc=para1_pairInf.get(analysed_pair_r);
		    					}
		    					neighbour_pairs=para1_bloclist.get(analysed_bloc);
			    				neighbour_pairs_list=neighbour_pairs.split(sep);
			    				nbr_synt_bloc=neighbour_pairs_list.length/2;
			    				med_ks=0.0;
		    					if (neighbour_pairs_list.length>0){
				    				int err_count=0;
				    				for (String p : neighbour_pairs_list){
				    					try{
				    						med_ks=med_ks+para1List.get(p);
				    					}catch (NullPointerException e){
				    						err_count=err_count+1;
				    					}
				    				}
				    				med_ks=med_ks/(neighbour_pairs_list.length);
			    				}else {
			    					med_ks=0.0;
			    				}
		    					
		    				}catch (NullPointerException e){
		    					KS=-1;
		    					nbr_synt_bloc=0;
		    					med_ks=0.0;
		    					System.out.println("nope");
		    				}
		    				resultFile.write(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
		    				resultFile.flush();
		    				//System.out.println(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
		    			}else if(para2_pairInf.containsKey(analysed_pair) || para2_pairInf.containsKey(analysed_pair_r)){
		    				type="inParalog";
		    				try {
		    					if(para2_pairInf.containsKey(analysed_pair)){
		    						try{
		    							KS=para2List.get(analysed_pair); // recuperer le Ks de la paire
				    				}catch(NullPointerException e){
				    					KS=-1; // recuperer le Ks de la paire
				    				}
				    				analysed_bloc=para2_pairInf.get(analysed_pair);
		    					}else{
		    						KS=para2List.get(analysed_pair_r);
				    				analysed_bloc=para2_pairInf.get(analysed_pair_r);
		    					}
			    				neighbour_pairs=para2_bloclist.get(analysed_bloc);
			    				neighbour_pairs_list=neighbour_pairs.split(sep);
			    				nbr_synt_bloc=neighbour_pairs_list.length/2;
			    				med_ks=0.0;
			    				if (neighbour_pairs_list.length>0){
				    				int err_count=0;
				    				for (String p : neighbour_pairs_list){
				    					try{
				    						med_ks=med_ks+para2List.get(p);
				    					}catch (NullPointerException e){
				    						err_count=err_count+1;
				    					}
				    				}
				    				med_ks=med_ks/(neighbour_pairs_list.length);
			    				}else {
			    					med_ks=0.0;
			    				}
		    				}catch (NullPointerException e){
		    					KS=-1;
		    					nbr_synt_bloc=0;
		    					med_ks=0.0;
		    					System.out.println("nope");
		    				}
		    				if(type!="")
							{
		    					resultFile.write(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
							    resultFile.flush();
							}
		    				//System.out.println(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n");
						}else if(tandem1_list.contains(analysed_pair) || tandem1_list.contains(analysed_pair_r)){
		    				type="tandem";		    				
							if(para2List.containsKey(analysed_pair)){
								KS=para2List.get(analysed_pair);
								analysed_bloc=para2_pairInf.get(analysed_pair);
							}
							else if (para2List.containsKey(analysed_pair_r)){
								KS=para2List.get(analysed_pair_r);
								analysed_bloc=para2_pairInf.get(analysed_pair_r);
							}
							else if (para1List.containsKey(analysed_pair)){
								KS=para1List.get(analysed_pair);
								analysed_bloc=para1_pairInf.get(analysed_pair);
							}
							else if (para1List.containsKey(analysed_pair_r)){
								KS=para1List.get(analysed_pair_r);
								analysed_bloc=para1_pairInf.get(analysed_pair_r);
							}
		    				nbr_synt_bloc=0;
		    				med_ks=0.0;
		    				if(type!="")
							{
		    					resultFile.write(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n"); 
							    resultFile.flush();
							}
						}else if(tandem2_list.contains(analysed_pair) || tandem2_list.contains(analysed_pair_r)){
		    				type="tandem";		    				
							if(para2List.containsKey(analysed_pair)){
								KS=para2List.get(analysed_pair);
								analysed_bloc=para2_pairInf.get(analysed_pair);
							}
							else if (para2List.containsKey(analysed_pair_r)){
								KS=para2List.get(analysed_pair_r);
								analysed_bloc=para2_pairInf.get(analysed_pair_r);
							}
							else if (para1List.containsKey(analysed_pair)){
								KS=para1List.get(analysed_pair);
								analysed_bloc=para1_pairInf.get(analysed_pair);
							}
							else if (para1List.containsKey(analysed_pair_r)){
								KS=para1List.get(analysed_pair_r);
								analysed_bloc=para1_pairInf.get(analysed_pair_r);
							}
		    				nbr_synt_bloc=0;
		    				med_ks=0.0;
		    				if(type!="")
							{
		    					resultFile.write(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n"); 
							    resultFile.flush();
							}
						}else if(pairList.containsKey(analysed_pair) || pairList.containsKey(analysed_pair_r)){
		    				type="outPara/miscellaneous";
		    				if(pairList.containsKey(analysed_pair)){
		    					KS=pairList.get(analysed_pair);
			    				analysed_bloc=pairInf.get(analysed_pair);
	    					}else{
	    						KS=pairList.get(analysed_pair_r);
	    	    				analysed_bloc=pairInf.get(analysed_pair_r);
	    					}
		    				neighbour_pairs=bloclist.get(analysed_bloc);
		    				neighbour_pairs_list=neighbour_pairs.split(sep);
		    				nbr_synt_bloc=neighbour_pairs_list.length/2;
		    				med_ks=0.0;
		    				if (neighbour_pairs_list.length>0){
			    				int err_count=0;
			    				for (String p : neighbour_pairs_list){
			    					try{
			    						med_ks=med_ks+pairList.get(p);
			    					}catch (NullPointerException e){
			    						err_count=err_count+1;
			    					}
			    				}
			    				med_ks=med_ks/(neighbour_pairs_list.length);
		    				}else {
		    					med_ks=0.0;
		    				}
		    				if(type!="")
							{
		    					resultFile.write(family_list.get(fam_count)+"\t"+family_list.get(t)+"\t"+type+"\t"+KS+"\t"+med_ks+"\t"+nbr_synt_bloc+"\n"); 
							    resultFile.flush();
							}
						}else{
		    				type="un";
		    				KS=-1;
		    				nbr_synt_bloc=0;
		    				med_ks=0.0;
		    			}

	    			}


	    		}
	    		fam_count ++;
	    	}

			
	    }catch(FileNotFoundException e){
	    	// fichier non trouve
	    	e.printStackTrace();
	    }catch (IOException e) {
	        // erreur d ecriture ou de lecture
	        e.printStackTrace();
	    }finally {
	    	
	        try {
	        	if (paires_ks != null)
	        		paires_ks.close();
	        	if (family != null)
	        		family.close();
	        	if (orthologous != null)
	        		orthologous.close();
	        	if (inparalogous1 != null)
	        		inparalogous1.close();
	        	if (inparalogous2 != null)
	        		inparalogous2.close();
	        	if (resultFile != null)
	        		resultFile.close();
	        } catch (IOException e) {
	           e.printStackTrace();
	        }
	    }
	    
	}
	  

	public static ArrayList<String> gettandem(String name){
		// methode Ks
//		
		String gene1=null;
    	String gene2=null;
    	String[] temp_array1;
    	String[] temp_array2;
    	String sep = "[|]+";
    	BufferedReader paires = null;
    	ArrayList<String> pairList = new ArrayList<String>(); //contient paires
    	
    	try{
    		paires = new BufferedReader(new FileReader(name));
        	String line;
        	String pair1= new String();
        	String pair2= new String();
        	
        	
        	while ((line = paires.readLine()) != null) {
        		String[] tand_array=line.split("\t");
        		pair1= new String();
        		pair2= new String();
        		if (tand_array.length >= 2){
        			gene1=tand_array[0].trim();
        			gene1=gene1.split(sep)[3];
        			gene1=gene1.toUpperCase();
        		    gene2=tand_array[1].trim();
        		    gene2=gene2.split(sep)[3];
        		    gene2=gene2.toUpperCase();
					
        			pair1=gene1+","+gene2;
        			pair2=gene2+","+gene1;
        			//System.out.println(gene1);
        			//System.out.println(gene2);
    				pairList.add(pair1);
    				pairList.add(pair2);
    				//System.out.println(pair1);
    				//System.out.println(pair2);

        		}
        		
        	}
        	return pairList;
        	
    	}catch(FileNotFoundException e){
	    	// fichier non trouve
	    	e.printStackTrace();
	    }catch (IOException e) {
	        // erreur d ecriture ou de lecture
	        e.printStackTrace();
	    }finally {
	        try {
	        	if (paires != null)
	        		paires.close();
	        	
	        } catch (IOException e) {
	           e.printStackTrace();
	        }
	    }
		return pairList;
	}
	
	
	public static HashMap<String, Double> getks(String name){
		// methode Ks
//		
		String geneks1=null;
    	String geneks2=null;
    	Double Ks;
    	String[] temp_array1;
    	String[] temp_array2;
    	String sep = "[|]+";
    	BufferedReader paires_ks = null;
    	HashMap<String, Double> pairList = new HashMap<String, Double>(); //contient paires et Ks
    	System.out.println("Find KS");
    	try{
    		paires_ks = new BufferedReader(new FileReader(name));
        	String line;
        	String pair1= new String();
        	String pair2= new String();
        	
        	
        	while ((line = paires_ks.readLine()) != null) {
        		String[] ks_array=line.split(",");
        		pair1= new String();
        		pair2= new String();
        		if (ks_array.length >= 2 && ks_array[0].contains(";")){
        			 geneks1=ks_array[0].split(";")[0].trim();
        			 // System.out.println(geneks1);
        			 geneks1=geneks1.split(sep)[3];
        			 geneks1=geneks1.toUpperCase();
        		     geneks2=ks_array[0].split(";")[1].trim();
        		     geneks2=geneks2.split(sep)[3];
        		     geneks2=geneks2.toUpperCase();
        		     Ks=Double.parseDouble(ks_array[1].trim());
        		     //System.out.println(geneks1);
        		     // temp_array1=geneks1.split("_");
 				    // if (temp_array1.length>=2){
				    	// System.out.println(temp_array[temp_array.length-1].split("."));
 				    	// temp_array1[temp_array1.length-1]=temp_array1[temp_array1.length-1].split("\\.")[0];
 				    	// temp_array1[temp_array1.length-1]=temp_array1[temp_array1.length-1].replace("T", "P");
 				    	// geneks1 = String.join("_", temp_array1);
 				    // }
 				    
 				   // temp_array2=geneks2.split("_");
				    // if (temp_array2.length>=2){
				    	// System.out.println(temp_array[temp_array.length-1].split("."));
				    	// temp_array2[temp_array2.length-1]=temp_array2[temp_array2.length-1].split("\\.")[0];
				    	// temp_array2[temp_array2.length-1]=temp_array2[temp_array2.length-1].replace("T", "P");
				    	// geneks2 = String.join("_", temp_array2);
				    // }

        			pair1=geneks1+","+geneks2;
        			pair2=geneks2+","+geneks1;
        			//System.out.println(pair1+" "+pair2);	
    				pairList.put(pair1, Ks);
    				pairList.put(pair2, Ks);
        		}
        		
        	}
        	return pairList;
        	
    	}catch(FileNotFoundException e){
	    	// fichier non trouve
	    	e.printStackTrace();
	    	System.out.println("File not found");
	    }catch (IOException e) {
	        // erreur d ecriture ou de lecture
	        e.printStackTrace();
	        System.out.println("Error");
	    }finally {
	        try {
	        	if (paires_ks != null)
	        		paires_ks.close();
	        	
	        } catch (IOException e) {
	           e.printStackTrace();
	        }
	    }
		return pairList;
	}
	
	public static ArrayList<HashMap<String, String>> parseblocfile(String name){
		BufferedReader inparalogous1_bloc = null;
		ArrayList<HashMap<String, String>> finallist= new ArrayList<HashMap<String, String>>();
		try {
			inparalogous1_bloc  = new BufferedReader(new FileReader(name));
			String line_para1;
//	    	double ks1;
	    	
	    	HashMap<String, String> genelist =new HashMap<String, String>(); // contient  nom et info genes
	    	HashMap<String, String> pairInf = new HashMap<String, String>(); //contient paires et bloc
	    	HashMap<String, String> syntlist =new HashMap<String, String>(); //contient un gene, et liste gene synt
	    	HashMap<String, String> bloclist =new HashMap<String, String>(); //contient un bloc, et liste genes compris

	    	String bloc = new String();
	    	String[] gene_array1 = new String[9];
	    	String[] gene_array2 = new String[9];
	    	String[] temp_array1 ;
	    	String[] temp_array2 ;
	    	String gene1_info = new String();
	    	String gene2_info = new String();
	    	String pair1 = new String();
	    	String pair2 = new String();
	    	String newvalue= new String();
	    	String cont_bloc= new String();
//	    	String line=new String();
	    	String chr1=new String();
//	    	String infgene1 =new String();
	    	String chr2=new String();
//	    	String infgene2 =new String();
	    	String gene_name1=new String();
	    	String gene_name2=new String();
	    	String sep = "[|]+";
	    	bloc = new String();
	    	
	   	
	    	while ((line_para1 = inparalogous1_bloc.readLine()) != null) {
	    		gene_array1 = new String[9];
	    		gene_array2 = new String[9];
	    		gene1_info = new String();
	    		gene2_info = new String();
	    		pair1 = new String();
	    		pair2 = new String();
	    		newvalue= new String();
				cont_bloc= new String();
//				line=new String();
				chr1=new String();
//				infgene1 =new String();
				chr2=new String();
//				infgene2 =new String();
				gene_name1=new String();
				gene_name2=new String();
				
	    		String[] line_array=line_para1.split("\t");
	    		
	    		if (line_array.length == 6 && line_array[0].startsWith("#")) // recup du num de bloc syntenique
	    		{
	    			String chrom1= line_array[2].trim();
	    			String chrom2= line_array[3].trim();
	    			String numbloc=line_array[0].substring(1);
	    			bloc = new String();
	    			bloc= numbloc+"_"+chrom1+"_"+chrom2;
	    			//System.out.println(bloc);
	    		}
	    		
	    		if (line_array.length >= 9 && !line_array[0].startsWith("#") ) //recup infos
		    	{
	    			chr1=line_array[0].trim();
	    			gene1_info=line_array[1].trim();
	    			gene1_info= gene1_info+"||"+bloc ;
	    			chr2=line_array[4].trim();
	    			gene2_info=line_array[5].trim();
	    			gene2_info= gene2_info+"||"+bloc ;
	    			
	    			gene_array1=gene1_info.split(sep);
	    			
	    			gene_name1=gene_array1[3].trim();
	    			gene_name1=gene_name1.toUpperCase();
					// temp_array1=gene_name1.split("_");
					
					// if (temp_array1.length>=2){
						// temp_array1[temp_array1.length-1]=temp_array1[temp_array1.length-1].split("\\.")[0];
						// temp_array1[temp_array1.length-1]=temp_array1[temp_array1.length-1].replace("T", "P");
						// gene_name1 = String.join("_", temp_array1);
					// }
	    			
	    			
	    			
	    			gene_array2=gene2_info.split(sep);
	    			gene_name2=gene_array2[3].trim();
	    			gene_name2=gene_name2.toUpperCase();
	    			
	    			// temp_array2=gene_name2.split("_");
					
					// if (temp_array2.length>=2){
						// temp_array2[temp_array2.length-1]=temp_array2[temp_array2.length-1].split("\\.")[0];
						// temp_array2[temp_array2.length-1]=temp_array2[temp_array2.length-1].replace("T", "P");
						// gene_name2 = String.join("_", temp_array2);
					// }
	    			
	    			pair1= gene_name1+","+gene_name2;
	    			pair2= gene_name2+","+gene_name1;
	    			
	    			//--------------------------------------------------------------Genes et infos -----------------------------------------------------------//
	    			
	    			if(genelist.containsKey(gene_name1))
	    			{

	    			}else
	    			{
	    				genelist.put(gene_name1, gene1_info);
	    			}
	    			if(genelist.containsKey(gene_name2))
	    			{

	    			}else
	    			{
	    				genelist.put(gene_name2, gene2_info);

	    			}
	    			
	    			
	    			
	    			//--------------------------------------------------------------Listes de genes synteniques a un gene -----------------------------------------------------------//

	    			if(syntlist.containsKey(gene_name1))
	    			{
	    				newvalue = syntlist.get(gene_name1);
	    				if(!newvalue.contains(gene_name2)){
	    					newvalue = newvalue+","+gene_name2;
		    				syntlist.put(gene_name1, newvalue);
	    				}

	    			}else
	    			{
	    				newvalue = new String();
	    				newvalue= gene_name2;
	    				syntlist.put(gene_name1, newvalue);

	    			}
	    			
	    			
	    			if(syntlist.containsKey(gene_name2))
	    			{
	    				newvalue = syntlist.get(gene_name2);
	    				if(!newvalue.contains(gene_name1)){
	    					newvalue = newvalue+","+gene_name1;
		    				syntlist.put(gene_name2, newvalue);
	    				}
	    			}else
	    			{
	    				newvalue= new String();
	    				newvalue=gene_name1;
	    				syntlist.put(gene_name2, newvalue);	    				
	    			}
	    			

	    			//--------------------------------------------------------------Contenu des blocs Synteniques -----------------------------------------------------------//
	    			
	    			if(bloclist.containsKey(bloc))
	    			{
	    				cont_bloc = bloclist.get(bloc);
	    				cont_bloc = cont_bloc+"||"+pair1;
	    				cont_bloc = cont_bloc+"||"+pair2;
	    				bloclist.put(bloc, cont_bloc);
	    			}else
	    			{
	    				cont_bloc = new String();
	    				cont_bloc= cont_bloc+"||"+pair1;
	    				cont_bloc= cont_bloc+"||"+pair2;
	    				bloclist.put(bloc, cont_bloc);
	    			}
	    			
	    			
	    			//--------------------------------------------------------------paires et bloc -----------------------------------------------------------//
	    			
	    			if(!pairInf.containsKey(pair1)){
	    				pairInf.put(pair1, bloc);
	    			}
	    			if(!pairInf.containsKey(pair2)){
	    				pairInf.put(pair2, bloc);
	    			}
//	   
	    		}
	    	}

	    	finallist.add(genelist);
	    	finallist.add(pairInf);
	    	finallist.add(syntlist);
	    	finallist.add(bloclist);
	    	
	    	
		}catch(FileNotFoundException e){
	    	// fichier non trouve
	    	e.printStackTrace();
	    }catch (IOException e) {
	        // erreur d ecriture ou de lecture
	        e.printStackTrace();
	    }finally {
	        try {
	        	if (inparalogous1_bloc != null)
	        		inparalogous1_bloc.close();
	        	
	        } catch (IOException e) {
	           e.printStackTrace();
	        }
	    }
		return finallist;
			
    	}
}


