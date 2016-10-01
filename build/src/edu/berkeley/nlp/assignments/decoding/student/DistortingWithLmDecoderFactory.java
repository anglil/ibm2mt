package edu.berkeley.nlp.assignments.decoding.student;


import java.util.*;
//import java.lang.reflect.Array;

import edu.berkeley.nlp.langmodel.NgramLanguageModel;
import edu.berkeley.nlp.mt.decoder.Decoder;
import edu.berkeley.nlp.mt.decoder.DecoderFactory;
import edu.berkeley.nlp.mt.decoder.DistortionModel;
import edu.berkeley.nlp.langmodel.EnglishWordIndexer;
//import edu.berkeley.nlp.mt.decoder.MonotonicGreedyDecoder;
//import edu.berkeley.nlp.mt.decoder.assert;
import edu.berkeley.nlp.mt.phrasetable.EnglishPhrase;
import edu.berkeley.nlp.mt.phrasetable.PhraseTable;
import edu.berkeley.nlp.mt.phrasetable.PhraseTableForSentence;
import edu.berkeley.nlp.mt.phrasetable.ScoredPhrasePairForSentence;
import edu.berkeley.nlp.util.GeneralPriorityQueue;
//import edu.berkeley.nlp.util.CollectionUtils;
import edu.berkeley.nlp.util.StringIndexer;
//import sun.reflect.generics.tree.Tree;





public class DistortingWithLmDecoderFactory implements DecoderFactory
{
	// implement the newDecoder interface defined in DecoderFactory
	public Decoder newDecoder(PhraseTable tm, NgramLanguageModel lm, DistortionModel dm) {
		
		// the distorting with lm decoder class that implements Decoder
		class DistortingWithLmDecoder implements Decoder
		{
			private PhraseTable tm;

			private NgramLanguageModel lm;

			private DistortionModel dm;
			
			/**
			 * 
			 * @param tm
			 * @param lm
			 * @param dm
			 */
			public DistortingWithLmDecoder(PhraseTable tm, NgramLanguageModel lm, DistortionModel dm){
				super();
				this.tm = tm;
				this.lm = lm;
				this.dm = dm;
			}
			
			class translationState
			{
				// state members
				private String lastWord1;
				private String lastWord2;
				private boolean[] bitString;
				private int endPointOfLastPhrase;
				
				// a score does not distinguish two states
				private double score;
				
				// back pointers
				private ScoredPhrasePairForSentence backPointerPhrase;
				private translationState backPointerState;
				
				
				
				/**
				 * 
				 * @param lastWord1
				 * @param lastWord2
				 * @param startByte
				 * @param endPointOfLastPhrase
				 * @param score
				 */
				public translationState(String lastWord1, String lastWord2,
						boolean[] startByte, int endPointOfLastPhrase, double score) {
					this.lastWord1 = lastWord1;
					this.lastWord2 = lastWord2;
					int bitStringLen = startByte.length;
					this.bitString = new boolean[bitStringLen]; // pass the reference
					System.arraycopy(startByte, 0, this.bitString, 0, bitStringLen);
					this.endPointOfLastPhrase = endPointOfLastPhrase;
					this.score = score;
					this.backPointerPhrase = null;
					this.backPointerState = null;
				}
				
				
				

				
				/**
				 * get the score of the state
				 * @return the score of the state
				 */
				public double getScore() {
					return score;
				}
				
				/**
				 * the unique string for every unique state
				 * @return
				 */
				public String getHashCodeVal() {
					String hashCodeVal = lastWord1 + lastWord2 + Arrays.toString(bitString) + Integer.toString(endPointOfLastPhrase);
					return hashCodeVal;
				}

//				/**
//				 * get the first one of the last two words in the translation 
//				 * @return the first one of the last two words in the translation
//				 */
//				public String getLastWord1() {
//					return lastWord1;
//				}
//				
//				/**
//				 * get the second one of the last two words in the translation 
//				 * @return the second one of the last two words in the translation
//				 */
//				public String getLastWord2() {
//					return lastWord2;
//				}
//				
//				/**
//				 * get the bit string in the translation state
//				 * @return
//				 */
//				public boolean[] getBitString() {
//					return bitString;
//				}
//				
//				/**
//				 * get the end point of the last phrase in the translation
//				 * @return the end point of the last phrase in the translation
//				 */
//				public double getEndPointOfLastPhrase() {
//					return endPointOfLastPhrase;
//				}
				
				/**
				 * get the number of 1's in the bit string
				 * @return the number of 1's in the bit string
				 */
				public int getBitTranslatedLen() {
					int count = 0;
					for (int bitStringWalker = 0; bitStringWalker < bitString.length; bitStringWalker++){
						if (bitString[bitStringWalker]) {
							count ++;
						}
					}
					return count;
				}

				

				/**
				 * get the phrase in the back pointer of this state
				 * @return
				 */
				public ScoredPhrasePairForSentence getBackPointerPhrase() {			
					return backPointerPhrase;
				}

				/**
				 * get the state in the back pointer of this state
				 * @return
				 */
				public translationState getBackPointerState() {
					return backPointerState;
				}
				
				/**
				 * 
				 * @param state
				 * @param phrase
				 */
				public void setBackPointer(translationState state, ScoredPhrasePairForSentence phrase) {
					backPointerPhrase = phrase;
					backPointerState = state;
				}

				/** 
				 * the state formed by combining state q with phrase p
				 * @param p a phrase to be combined with this state
				 * @return a new state
				 */
				public translationState nextState(ScoredPhrasePairForSentence p, PhraseTableForSentence phraseTable) {
					
					String epsilon1Neg = new String(lastWord1);
					String epsilon0 = new String(lastWord2);
					StringIndexer engInd = EnglishWordIndexer.getIndexer();
					int e1 = engInd.addAndGetIndex(epsilon1Neg);
					int e2 = engInd.addAndGetIndex(epsilon0);
					
					
					List<String> engWords = p.getEnglish();
					int engLen = engWords.size(); 
					String e1Prime, e2Prime;
					if (engLen >= 2) {
						e1Prime = new String(engWords.get(engLen-2));
						e2Prime = new String(engWords.get(engLen-1));
					} else {
						e1Prime = new String(lastWord1);
						e2Prime = new String(engWords.get(engLen-1));
					}
				
					int s = p.getStart();
					int t = p.getEnd();
					int foreignLen = bitString.length;
					boolean[] bitStringPrime = new boolean[foreignLen];
					System.arraycopy(bitString, 0, bitStringPrime, 0, foreignLen);
					
					// update the bit string
					for (int foreignSentWalker = 0; foreignSentWalker < foreignLen; foreignSentWalker++){
						if ((foreignSentWalker >= s-1) && (foreignSentWalker <= t-1)) {
							bitStringPrime[foreignSentWalker] = true;
						}
					}
					int endPointOfLastPhrasePrime = t;
					
					// attach the first two words to the language model
					int[] indexedEnglishSent = p.english.indexedEnglish;
					int lenSent = indexedEnglishSent.length;
					
					int[] extendedEnglishSent = new int[lenSent+2];
					extendedEnglishSent[0] = e1;
					extendedEnglishSent[1] = e2;
					for (int sentWalker = 2; sentWalker < lenSent+2; sentWalker++) {
						extendedEnglishSent[sentWalker] = indexedEnglishSent[sentWalker-2];
					}
					
					double scorePrime = score + scoreTmLmDm(lm.getOrder(), p.getForeignLength(), extendedEnglishSent, 
							lenSent, lm, phraseTable, dm, s, t, endPointOfLastPhrase, s, p); 
					
					return new translationState(e1Prime, e2Prime, bitStringPrime, endPointOfLastPhrasePrime, scorePrime);
				}
				
				
//				/**
//				 * compare two states by comparing the first four parameters in the state
//				 * @param qPrime
//				 * @return
//				 */
//				public Boolean eq(translationState qPrime) {
//					if (qPrime.getLastWord1().equals(lastWord1) && qPrime.getLastWord2().equals(lastWord2) && qPrime.getBitString().equals(bitString) && qPrime.getEndPointOfLastPhrase() == endPointOfLastPhrase) {
//						return true;
//					}
//					return false;
//				}
				
				
				/**
				 * select phrases
				 * @return
				 */
				public List<ScoredPhrasePairForSentence> ph(List<ScoredPhrasePairForSentence> allPhrases) {
					
					List<ScoredPhrasePairForSentence> selectPhrases = new ArrayList<ScoredPhrasePairForSentence>();
					//PhraseTableForSentence tmState = tm.initialize(sentence);
					//List<ScoredPhrasePairForSentence> phrases = tmState.getScoreSortedTranslationsForSpan(0, tmState.getMaxPhraseLength());
					for (final ScoredPhrasePairForSentence phrase : allPhrases) {
						int start = phrase.getStart();
						int end = phrase.getEnd();
						if (!ifOverlap(start, end) && !ifDistortionLimitViolated(phrase)) {
							selectPhrases.add(phrase);
						}
					}
					return selectPhrases;
				}
				
				
				/**
				 * 
				 * @param start
				 * @param end
				 * @return
				 */
				public Boolean ifOverlap(int start, int end) {
					for (int bitStringWalker = start; bitStringWalker <= end; bitStringWalker++) {
						if (bitString[bitStringWalker]) {
							return true;
						}
					}
					return false;
				}
				
				/**
				 * 
				 * @param p
				 * @return
				 */
				public Boolean ifDistortionLimitViolated(ScoredPhrasePairForSentence p) {
					if (Math.abs(endPointOfLastPhrase + 1 - p.getStart()) <= dm.getDistortionLimit()) {
						return false;
					}
					return true;
				}
				
			}
			// end of the state class
			
			
			/**
			 * input: tokenized foreign sentence
			 * output: phrase-based English sentence
			 */
			public List<ScoredPhrasePairForSentence> decode(List<String> sentence) {
				int length = sentence.size();
				// phrase table initialization
				PhraseTableForSentence tmState = tm.initialize(sentence);
				
				List<ScoredPhrasePairForSentence> allPhrases = new ArrayList<ScoredPhrasePairForSentence>();
				
				for (int i = 0; i<=sentence.size()-1; i++) {
					for (int j=1; j < Math.min(tmState.getMaxPhraseLength(), sentence.size()-i); j++) {
						List<ScoredPhrasePairForSentence> coveredPhrases = tmState.getScoreSortedTranslationsForSpan(i, i+j);
						if (coveredPhrases != null) {
							allPhrases.addAll(coveredPhrases);
						}
					}
				}
				
				// initialization of the set of states
				ArrayList<TreeMap<String, translationState>> stateSetAll = new ArrayList<TreeMap<String, translationState>>(length+1);
				ArrayList<GeneralPriorityQueue<translationState>> stateSetAllBeam = new ArrayList<GeneralPriorityQueue<translationState>>(length+1);
				for (int stateSetWalker = 0; stateSetWalker <= length; stateSetWalker++) {
					stateSetAll.add(new TreeMap<String, translationState>());
					stateSetAllBeam.add(new GeneralPriorityQueue<translationState>());
//					stateSetAll.set(stateSetWalker, new TreeMap<String, translationState>());
				}
				
				
				// the first state in the first state set
				boolean[] startByte = new boolean[length]; 
				for (int byteWalker = 0; byteWalker < length; byteWalker++) {
					startByte[byteWalker] = false;
				}
				translationState state0 = new translationState(NgramLanguageModel.START, NgramLanguageModel.START, startByte, 0, 0.0);

				
				// the state set that includes the hash code of a state
				// the hash map part
				TreeMap<String, translationState> stateSet0 = new TreeMap<String, translationState>();
				stateSet0.put(state0.getHashCodeVal(), state0);
				stateSetAll.set(0, stateSet0);
				// the priority queue part
				GeneralPriorityQueue<translationState> stateSetFirstBeam = new GeneralPriorityQueue<translationState>();
				stateSetFirstBeam.setPriority(state0, -state0.getScore()); // the priority is reversed
				stateSetAllBeam.set(0, stateSetFirstBeam);
				
				
				
//				System.out.println("Beginning Decoding");
				
				
				// the decoding algorithm
				for (int stateSetWalker = 0; stateSetWalker < length; stateSetWalker++) {
					// print the index of the word the translator is arriving at
					System.out.println(stateSetWalker);
					
					
					// compute the beam of a state set
//					TreeMap<String, translationState> stateSetBeam = beam(stateSetAll.get(stateSetWalker));
					TreeMap<String, translationState> stateSetBeam = stateSetAll.get(stateSetWalker); // the beam sized is kept within 2000, so no need for a beam search
					
					
//					System.out.println("stateSetBeam size:");
//					System.out.println(stateSetBeam.size());
//					int ctrOut = 0;
					for (final Map.Entry<String, translationState> kvp: stateSetBeam.entrySet()) {
//						System.out.println("Enter New State");
//						System.out.println(ctrOut);
//						ctrOut++;
//						System.out.println(stateSetBeam.size());
						// get the phrase set of a translation state
						//int ctrIn = 0;
						translationState q = kvp.getValue();
						List<ScoredPhrasePairForSentence> phraseSet = q.ph(allPhrases);
						
						for (final ScoredPhrasePairForSentence p : phraseSet) {
							//System.out.println("Entered Add Phrase");
							//System.out.println(ctrIn);
							//System.out.println(phraseSet.size());
							//ctrIn++;
							// get the next translation state
							
							translationState qPrime = q.nextState(p, tmState);
							
							// update the state set of the j-th element of the state set list
//							int j = qPrime.getBitTranslatedLen();
							
//							List<translationState> stateSetJ = stateSetAll.get(j);
//							stateSetAll.set(j, addState(stateSetJ, qPrime, q, p));
							
							// the number of ones in the bit string
							int j = qPrime.getBitTranslatedLen(); 
							addStateToStateSet(stateSetAll.get(j), stateSetAllBeam.get(j), qPrime, q, p);
							
							//System.out.println(ctrIn);
							//System.out.println(phraseSet.size());
							//System.out.println("Exit Add Phrase");
							
						}
//						System.out.println("Exit New State");
//						System.out.println(ctrOut);
//						System.out.println(stateSetBeam.size());
						
						
					}
				}
				
				System.out.println("Completed Decoding");
				
				
				// get the state with the highest score
				int stateSetSize = length;
				boolean flag = true;
				GeneralPriorityQueue<translationState> stateSetNBeam = new GeneralPriorityQueue<translationState>();
				while (flag) { 
					stateSetNBeam = stateSetAllBeam.get(stateSetSize-1);
					if (stateSetNBeam.size() != 0) {
						
						flag = false;
					}
					stateSetSize--; 
				}
//				GeneralPriorityQueue<translationState> stateSetNBeam = stateSetAllBeam.get(length); // the last state set in stateSetAll
				
				
				// get the highest scoring state
				translationState stateHighestScore = null;
				while (!stateSetNBeam.isEmpty()) {
					stateHighestScore = stateSetNBeam.removeFirst();
				}
				
//				translationState stateHighestScore = null;
////				double maxScore = Double.NEGATIVE_INFINITY;
//				for (final Map.Entry<String, translationState> kvp : stateSetN.entrySet()) {
////					double scoreTmp = stateTmp.getScore();
////					if (scoreTmp > maxScore) {
////						maxScore = scoreTmp;
////						stateHighestScore = stateTmp;
////					}
//					
//					translationState stateTmp = kvp.getValue();
//					stateHighestScore = stateTmp;
//				}
				
				System.out.println("Completed Best State Value Estimation");
				
				
				// from the state with the highest score, recover the highest scoring translation, by tracing the back-pointers for the states
				Stack<ScoredPhrasePairForSentence> res = new Stack<ScoredPhrasePairForSentence>();
				
//				List<ScoredPhrasePairForSentence> resSent = new ArrayList<ScoredPhrasePairForSentence>();
				while (stateHighestScore != null) {
					ScoredPhrasePairForSentence pHighestScore = stateHighestScore.getBackPointerPhrase();
					res.push(pHighestScore);
					stateHighestScore = stateHighestScore.getBackPointerState();
				}
				
				List<ScoredPhrasePairForSentence> resSent = new ArrayList<ScoredPhrasePairForSentence>();
				// pop out the null
				res.pop();
				while (!res.empty()) {
					resSent.add(res.pop());
				}
				
				int f = resSent.get(resSent.size()-1).getEnd();
				
				List<String> fs = new ArrayList<String>();
				// add a period to the end of the translated sentence
				String[] endPeriod = new String[1];
				endPeriod[0] = ".";
				fs.add(".");
				
				resSent.add(new ScoredPhrasePairForSentence(new EnglishPhrase(endPeriod), 1.0, fs, f, f+1));
				
				return resSent;
				
				
			}
			
			
			
			/**
			 * 
			 * @param stateSetJ
			 * @param qPrime
			 * @param q
			 * @param p
			 * @return 
			 * @return
			 */
			private void addStateToStateSet(TreeMap<String, translationState> stateSetJ, GeneralPriorityQueue<translationState> stateSetJBeam, translationState qPrime, translationState q, ScoredPhrasePairForSentence p) {
				boolean isStateExist = false;
				if (stateSetJ.size() != 0) {
					// check whether the new state is already in the current state set (eq)
					if (stateSetJ.containsKey(qPrime.getHashCodeVal())) {
						isStateExist = true;
						translationState qPrimePrime = stateSetJ.get(qPrime.getHashCodeVal());
						if (qPrime.getScore() > qPrimePrime.getScore()) {
							stateSetJ.remove(qPrimePrime.getHashCodeVal());
							stateSetJBeam.removeKey(qPrimePrime);
							qPrime.setBackPointer(q, p);
							stateSetJ.put(qPrime.getHashCodeVal(), qPrime);
							stateSetJBeam.setPriority(qPrime, -qPrime.getScore());
						}
					}
				}
				if (!isStateExist) {
					// keep beam size less than 2000
					if (stateSetJ.size() < 2000) {
						qPrime.setBackPointer(q, p);
						stateSetJ.put(qPrime.getHashCodeVal(), qPrime);
						stateSetJBeam.setPriority(qPrime, -qPrime.getScore());
					} else {
						translationState stateLowScore = stateSetJBeam.removeFirst(); // remove the one with the lowest score
						
						// if qPrime is not the lowest scoring one
						if (qPrime.getScore() > -stateLowScore.getScore()) {
							stateSetJ.remove(stateLowScore.getHashCodeVal());
							qPrime.setBackPointer(q, p);
							stateSetJ.put(qPrime.getHashCodeVal(), qPrime);
							stateSetJBeam.setPriority(qPrime, -qPrime.getScore());
						} else { // qPrime is the lowest scoring one
							stateSetJBeam.setPriority(stateLowScore, -stateLowScore.getScore());
						}
					}
				}
			}


			
			/**
			 * 
			 * @param lmOrder
			 * @param prevLmStateLength
			 * @param lmStateBuf
			 * @param totalTrgLength
			 * @param lm
			 * @param tm
			 * @param dm
			 * @param s
			 * @param t
			 * @param endPrev
			 * @param beginCurr
			 * @return
			 */
			private double scoreTmLmDm(final int lmOrder, final int prevLmStateLength, final int[] lmStateBuf, 
					final int totalTrgLength, final NgramLanguageModel lm, final PhraseTableForSentence tm, 
					final DistortionModel dm, final int s, final int t, final int endPrev, final int beginCurr, final ScoredPhrasePairForSentence phraseCurr) {
				double scoreLmUpdate = scoreLm(lmOrder, prevLmStateLength, lmStateBuf, totalTrgLength, lm);
//				List<ScoredPhrasePairForSentence> tmp = tm.getScoreSortedTranslationsForSpan(s, t);
//				double scoreTmUpdate = tmp.get(0).score;
				double scoreTmUpdate = phraseCurr.score;
				double scoreDmUpdate = dm.getDistortionScore(endPrev, beginCurr);
				return scoreLmUpdate + scoreTmUpdate + scoreDmUpdate;
			}
			
			
			/** language model score
			 * 
			 * @param lmOrder	language model order
			 * @param prevLmStateLength	previous language model state length
			 * @param lmStateBuf	language model state buffer
			 * @param totalTrgLength	total trigram length
			 * @param lm	language model
			 * @return
			 */
			private double scoreLm(final int lmOrder, final int prevLmStateLength, final int[] lmStateBuf, final int totalTrgLength, final NgramLanguageModel lm) {
				double score = 0.0;

				if (prevLmStateLength < lmOrder - 1) {
					for (int i = 1; prevLmStateLength + i < lmOrder; ++i) {
						final double lmProb = lm.getNgramLogProbability(lmStateBuf, 0, prevLmStateLength + i);
						score += lmProb;
					}
				}
				for (int i = 0; i <= totalTrgLength - lmOrder; ++i) {
					final double lmProb = lm.getNgramLogProbability(lmStateBuf, i, i + lmOrder);
					score += lmProb;
				}
				return score;
			}
			
			
		}
		
		System.out.println("distorting with lm decoder...");
		return new DistortingWithLmDecoder(tm, lm, dm);
	}

}





