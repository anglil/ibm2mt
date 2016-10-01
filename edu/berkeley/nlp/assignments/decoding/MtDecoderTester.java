package edu.berkeley.nlp.assignments.decoding;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.berkeley.nlp.assignments.decoding.student.MultiBeamDecoderFactory;
import edu.berkeley.nlp.assignments.decoding.student.DistortingWithLmDecoderFactory;
import edu.berkeley.nlp.assignments.decoding.student.MonotonicNoLmDecoderFactory;
import edu.berkeley.nlp.assignments.decoding.student.MonotonicWithLmDecoderFactory;
import edu.berkeley.nlp.io.IOUtils;
import edu.berkeley.nlp.langmodel.EnglishWordIndexer;
import edu.berkeley.nlp.langmodel.NgramLanguageModel;
import edu.berkeley.nlp.langmodel.RandomLanguageModel;
import edu.berkeley.nlp.langmodel.StubLanguageModel;
import edu.berkeley.nlp.langmodel.impl.KneserNeyLm;
import edu.berkeley.nlp.langmodel.impl.NgramLanguageModelAdaptor;
import edu.berkeley.nlp.langmodel.impl.NgramMapOpts;
import edu.berkeley.nlp.mt.BleuScore;
import edu.berkeley.nlp.mt.decoder.Decoder;
import edu.berkeley.nlp.mt.decoder.Decoder.StaticMethods;
import edu.berkeley.nlp.mt.decoder.DecoderFactory;
import edu.berkeley.nlp.mt.decoder.DistortionModel;
import edu.berkeley.nlp.mt.decoder.LinearDistortionModel;
import edu.berkeley.nlp.mt.decoder.Logger;
import edu.berkeley.nlp.mt.decoder.MonotonicGreedyDecoder.MonotonicGreedyDecoderFactory;
import edu.berkeley.nlp.mt.decoder.StubDistortionModel;
import edu.berkeley.nlp.mt.phrasetable.PhraseTable;
import edu.berkeley.nlp.mt.phrasetable.ScoredPhrasePairForSentence;
import edu.berkeley.nlp.util.CollectionUtils;
import edu.berkeley.nlp.util.CommandLineUtils;
import edu.berkeley.nlp.util.Counter;
import edu.berkeley.nlp.util.Pair;
import edu.berkeley.nlp.util.StrUtils;

/**
 * This is the main harness for assignment 2. To run this harness, use
 * <p/>
 * java edu.berkeley.nlp.assignments.assign1.MtDecoderTester -path
 * ASSIGNMENT_DATA_PATH -decoderType DECODER_TYPE_DESCRIPTOR_STRING
 * <p/>
 * 
 * 
 * @author Adam Pauls
 */
public class MtDecoderTester
{

	// HELPER CLASS FOR THE HARNESS, CAN IGNORE
	static class SentenceCollection implements Iterable<List<String>>
	{
		static class SentenceIterator implements Iterator<List<String>>
		{

			Iterator<String> reader;

			public boolean hasNext() {
				return reader.hasNext();
			}

			public List<String> next() {
				String line = reader.next();
				String[] words = line.split("\\s+");
				List<String> sentence = new ArrayList<String>();
				for (int i = 0; i < words.length; i++) {
					String word = words[i];
					sentence.add(word.toLowerCase());
				}
				return sentence;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}

			public SentenceIterator(Iterator<String> reader) {
				this.reader = reader;
			}
		}

		String fileName;

		public Iterator<List<String>> iterator() {
			try {

				return new SentenceIterator(IOUtils.lineIterator(fileName));
			} catch (FileNotFoundException e) {
				System.err.println("Problem with SentenceIterator for " + fileName);
				throw new RuntimeException(e);
			} catch (IOException e) {
				System.err.println("Problem with SentenceIterator for " + fileName);
				throw new RuntimeException(e);

			}
		}

		public SentenceCollection(String fileName) {
			this.fileName = fileName;
		}

		public static class Reader
		{
			static Iterable<List<String>> readSentenceCollection(String fileName) {
				return new SentenceCollection(fileName);
			}
		}

	}

	enum DecoderType
	{
		MONO_GREEDY
		{
			@Override
			public DecoderFactory getFactory() {
				return new MonotonicGreedyDecoderFactory();
			}

			@Override
			public DistortionModel getDistortionModel(Counter<String> weights) {
				return new StubDistortionModel();
			}

			@Override
			public NgramLanguageModel getLanguageModel(File lmFile, boolean random) {
				return getActualLanguageModel(lmFile, random);
			}
		},
		MONO_NOLM
		{
			@Override
			public DecoderFactory getFactory() {
				return new MonotonicNoLmDecoderFactory();
			}

			@Override
			public DistortionModel getDistortionModel(Counter<String> weights) {
				return new StubDistortionModel();
			}

			@Override
			public NgramLanguageModel getLanguageModel(File lmFile, boolean random) {
				return new StubLanguageModel();
			}
		},
		MONO_LM
		{
			@Override
			public DecoderFactory getFactory() {
				return new MonotonicWithLmDecoderFactory();
			}

			@Override
			public DistortionModel getDistortionModel(Counter<String> weights) {
				return new StubDistortionModel();
			}

			@Override
			public NgramLanguageModel getLanguageModel(File lmFile, boolean random) {
				return getActualLanguageModel(lmFile, random);
			}
		},
		DIST_LM
		{
			@Override
			public DecoderFactory getFactory() {
				return new DistortingWithLmDecoderFactory();
			}

			@Override
			public DistortionModel getDistortionModel(Counter<String> weights) {
				return getActualDistortionModel(weights);
			}

			@Override
			public NgramLanguageModel getLanguageModel(File lmFile, boolean random) {
				return getActualLanguageModel(lmFile, random);
			}

		},
		MULTIBEAM
		{
			@Override
			public DecoderFactory getFactory() {
				return new MultiBeamDecoderFactory();
			}

			@Override
			public DistortionModel getDistortionModel(Counter<String> weights) {
				return getActualDistortionModel(weights);
			}

			@Override
			public NgramLanguageModel getLanguageModel(File lmFile, boolean random) {
				return getActualLanguageModel(lmFile, random);
			}
		};

		public abstract DecoderFactory getFactory();

		public abstract DistortionModel getDistortionModel(Counter<String> weights);

		public abstract NgramLanguageModel getLanguageModel(File lmFile, boolean random);
	}

	public static void main(String[] args) {
		// Parse command line flags and arguments
		Map<String, String> argMap = CommandLineUtils.simpleCommandLineParser(args);

		// Set up default parameters and settings
		String basePath = ".";
		DecoderType decoderType = DecoderType.MONO_GREEDY;
		// You can use this to make decoding runs run in less time, but remember that we will
		// evaluate you on all test sentences.
		int maxNumTest = Integer.MAX_VALUE;
		boolean sanityCheck = false;
		boolean printTranslations = true;
		boolean randomLm = false;

		// Update defaults using command line specifications

		// The path to the assignment data
		if (argMap.containsKey("-path")) {
			basePath = argMap.get("-path");
		}
		System.out.println("Using base path: " + basePath);

		// A string descriptor of the model to use
		if (argMap.containsKey("-decoderType")) {
			decoderType = DecoderType.valueOf(argMap.get("-decoderType"));
		}
		System.out.println("Using decoderType: " + decoderType);

		if (argMap.containsKey("-maxNumTest")) {
			maxNumTest = Integer.parseInt(argMap.get("-maxNumTest"));
		}
		System.out.println("Decoding " + (maxNumTest == Integer.MAX_VALUE ? "all" : ("" + maxNumTest)) + " sentences.");

		if (argMap.containsKey("-noprint")) {
			printTranslations = false;
		}

		if (argMap.containsKey("-sanityCheck")) {
			sanityCheck = true;
		}
		if (sanityCheck) System.out.println("Only doing sanity check.");

		// Use an LM which just returns a random (but consistent for a given n-gram) score. This should speed up loading and facilitate sanity checking
		if (argMap.containsKey("-randomLm")) {
			randomLm = true;
		}

		if (sanityCheck) randomLm = true;
		String prefix = sanityCheck ? "sanity_" : "";

		// Read in all the assignment data
		File lmFile = new File(basePath, prefix + "lm.gz");
		File phraseTableFile = new File(basePath, prefix + "phrasetable.txt.gz");

		File testFrench = new File(basePath, prefix + "test.fr");
		File testEnglish = new File(basePath, prefix + "test.en");
		File weightsFile = new File(basePath, "weights.txt");

		// Build the language model
		DecoderFactory languageModelFactory = decoderType.getFactory();

		PhraseTable phraseTable = new PhraseTable(5, 30);
		final Counter<String> weights = readWeightsFile(weightsFile);
		final DistortionModel distortionModel = decoderType.getDistortionModel(weights);
		phraseTable.readFromFile(phraseTableFile.getPath(), weights);
		final NgramLanguageModel languageModel = decoderType.getLanguageModel(lmFile, randomLm);
		Decoder decoder = languageModelFactory.newDecoder(phraseTable, languageModel, distortionModel);
		evaluateDecoder(decoder, phraseTable, testFrench, testEnglish, weightsFile, getActualLanguageModel(lmFile, randomLm), maxNumTest, printTranslations,
			getActualDistortionModel(weights));

	}

	/**
	 * @param useTestSet
	 * @param phraseTableFile
	 * @param devFrench
	 * @param devEnglish
	 * @param testFrench
	 * @param testEnglish
	 * @param weightsFile
	 * @param languageModel
	 */
	private static void evaluateDecoder(Decoder decoder, PhraseTable phraseTable, File testFrench, File testEnglish, File weightsFile,
		NgramLanguageModel languageModel, int maxNumTest, boolean printTranslations, DistortionModel distortionModel) {

		printMemoryUsage();
		final String frenchData = (testFrench).getPath();
		Iterable<List<String>> frenchSentences = SentenceCollection.Reader.readSentenceCollection(frenchData);
		final String englishData = (testEnglish).getPath();
		Iterable<List<String>> englishSentences = SentenceCollection.Reader.readSentenceCollection(englishData);
		List<BleuScore> scores = new ArrayList<BleuScore>();
		double[] modelScore = new double[1];
		doDecoding(decoder, frenchSentences, englishSentences, scores, maxNumTest, printTranslations, languageModel, modelScore, distortionModel);
		String bleuString = new BleuScore(scores).toString();
		System.out.println("BLEU score on " + ("test") + " data was " + bleuString);
		System.out.println("Total model score on " + ("test") + " data was " + modelScore[0]);

	}

	/**
	 * 
	 */
	private static void printMemoryUsage() {
		System.gc();
		System.gc();
		System.gc();
		System.out.println("Memory usage is " + getUsedMemoryStr());
	}

	/**
	 * @param weights
	 * @return
	 */
	private static LinearDistortionModel getActualDistortionModel(Counter<String> weights) {
		return new LinearDistortionModel(4, weights.getCount("linearDist"));
	}

	private static NgramLanguageModel actualLanguageModel = null;

	/**
	 * @param lmFile
	 * @param random
	 * @return
	 */
	private static NgramLanguageModel getActualLanguageModel(File lmFile, boolean random) {
		if (actualLanguageModel == null)
			actualLanguageModel = random ? new RandomLanguageModel() : new NgramLanguageModelAdaptor(KneserNeyLm.fromFile(new NgramMapOpts(), lmFile.getPath(),
				3, EnglishWordIndexer.getIndexer()));
		return actualLanguageModel;
	}

	/**
	 * @param decoder
	 * @param frenchSentences
	 * @param englishSentences
	 * @param scores
	 * @param languageModel
	 */
	private static void doDecoding(Decoder decoder, Iterable<List<String>> frenchSentences, Iterable<List<String>> englishSentences, List<BleuScore> scores,
		int maxNumTest, boolean printTranslations, NgramLanguageModel languageModel, double[] modelScore, DistortionModel dm) {
		long startTime = System.nanoTime();
		int sent = 0;
		System.out.println("Decoding " + (maxNumTest == Integer.MAX_VALUE ? "all" : ("" + maxNumTest)) + " test sentences");
		for (Pair<List<String>, List<String>> input : zip(Pair.makePair(frenchSentences, englishSentences))) {
			if (sent >= maxNumTest) break;
			sent++;

			if (sent % 100 == 0) Logger.logs("On sentence " + sent);

			final List<ScoredPhrasePairForSentence> hyp = decoder.decode(input.getFirst());
			double score = StaticMethods.scoreHypothesis(hyp, languageModel, dm);
			List<String> hypothesisEnglish = Decoder.StaticMethods.extractEnglish(hyp);
			List<String> reference = input.getSecond();
			if (printTranslations) {

				System.out.println("Model score:\t" + score);
				System.out.println("Input:\t\t" + StrUtils.join(input.getFirst()));
				System.out.println("Hypothesis\t" + StrUtils.join(hypothesisEnglish));
				System.out.println("Reference:\t" + StrUtils.join(reference));
				System.out.println();
			}
			BleuScore bleuScore = new BleuScore(hypothesisEnglish, reference);
			modelScore[0] += score;
			scores.add(bleuScore);
		}
		long endTime = System.nanoTime();

		System.out.println("Decoding took " + BleuScore.formatDouble((endTime - startTime) / 1e9) + "s");
	}

	private static String bytesToString(long b) {
		double gb = (double) b / (1024 * 1024 * 1024);
		if (gb >= 1) return gb >= 10 ? (int) gb + "G" : round(gb, 1) + "G";
		double mb = (double) b / (1024 * 1024);
		if (mb >= 1) return mb >= 10 ? (int) mb + "M" : round(mb, 1) + "M";
		double kb = (double) b / (1024);
		if (kb >= 1) return kb >= 10 ? (int) kb + "K" : round(kb, 1) + "K";
		return b + "";
	}

	private static double round(double x, int numPlaces) {
		double scale = Math.pow(10, numPlaces);
		return Math.round(x * scale) / scale;
	}

	private static String getUsedMemoryStr() {
		long totalMem = Runtime.getRuntime().totalMemory();
		long freeMem = Runtime.getRuntime().freeMemory();
		return bytesToString(totalMem - freeMem);
	}

	private static <S, T> Iterable<Pair<S, T>> zip(Pair<? extends Iterable<S>, ? extends Iterable<T>> input) {
		final Iterator<S> firstIterator = input.getFirst().iterator();
		final Iterator<T> secondIterator = input.getSecond().iterator();
		return CollectionUtils.iterable(new Iterator<Pair<S, T>>()
		{

			public boolean hasNext() {
				return firstIterator.hasNext() && secondIterator.hasNext();
			}

			public Pair<S, T> next() {
				return Pair.makePair(firstIterator.next(), secondIterator.next());
			}

			public void remove() {

				throw new UnsupportedOperationException("Remove not supported");
			}
		});
	}

	private static Counter<String> readWeightsFile(File weightsFile) {
		Counter<String> ret = new Counter<String>();
		try {
			for (String line : CollectionUtils.iterable(IOUtils.lineIterator(weightsFile.getPath()))) {
				if (line.trim().length() == 0) continue;
				String[] parts = line.trim().split("\t");
				ret.setCount(parts[0].intern(), Double.parseDouble(parts[1]));

			}
		} catch (NumberFormatException e) {
			System.err.println("Error reading weights file " + weightsFile);
			throw new RuntimeException(e);

		} catch (IOException e) {
			System.err.println("Error reading weights file " + weightsFile);
			throw new RuntimeException(e);

		}
		return ret;
	}
}
