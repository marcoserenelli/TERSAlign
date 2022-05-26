// Generated from C:/Users/Filippo/Projects/TertiaryStructuresComparator/src/main/resources\TertiaryStructureBonds.g4 by ANTLR 4.10.1

	package it.unicam.cs.bdslab.tersaling;

import org.antlr.v4.runtime.atn.*;
import org.antlr.v4.runtime.dfa.DFA;
import org.antlr.v4.runtime.*;
import org.antlr.v4.runtime.misc.*;
import org.antlr.v4.runtime.tree.*;
import java.util.List;

@SuppressWarnings({"all", "warnings", "unchecked", "unused", "cast"})
public class TertiaryStructureBondsParser extends Parser {
	static { RuntimeMetaData.checkVersion("4.10.1", RuntimeMetaData.VERSION); }

	protected static final DFA[] _decisionToDFA;
	protected static final PredictionContextCache _sharedContextCache =
		new PredictionContextCache();
	public static final int
		T__0=1, T__1=2, T__2=3, T__3=4, INDEX=5, LINE_COMMENT=6, WS=7;
	public static final int
		RULE_bonds = 0, RULE_bond = 1;
	private static String[] makeRuleNames() {
		return new String[] {
			"bonds", "bond"
		};
	}
	public static final String[] ruleNames = makeRuleNames();

	private static String[] makeLiteralNames() {
		return new String[] {
			null, "';'", "'('", "','", "')'"
		};
	}
	private static final String[] _LITERAL_NAMES = makeLiteralNames();
	private static String[] makeSymbolicNames() {
		return new String[] {
			null, null, null, null, null, "INDEX", "LINE_COMMENT", "WS"
		};
	}
	private static final String[] _SYMBOLIC_NAMES = makeSymbolicNames();
	public static final Vocabulary VOCABULARY = new VocabularyImpl(_LITERAL_NAMES, _SYMBOLIC_NAMES);

	/**
	 * @deprecated Use {@link #VOCABULARY} instead.
	 */
	@Deprecated
	public static final String[] tokenNames;
	static {
		tokenNames = new String[_SYMBOLIC_NAMES.length];
		for (int i = 0; i < tokenNames.length; i++) {
			tokenNames[i] = VOCABULARY.getLiteralName(i);
			if (tokenNames[i] == null) {
				tokenNames[i] = VOCABULARY.getSymbolicName(i);
			}

			if (tokenNames[i] == null) {
				tokenNames[i] = "<INVALID>";
			}
		}
	}

	@Override
	@Deprecated
	public String[] getTokenNames() {
		return tokenNames;
	}

	@Override

	public Vocabulary getVocabulary() {
		return VOCABULARY;
	}

	@Override
	public String getGrammarFileName() { return "TertiaryStructureBonds.g4"; }

	@Override
	public String[] getRuleNames() { return ruleNames; }

	@Override
	public String getSerializedATN() { return _serializedATN; }

	@Override
	public ATN getATN() { return _ATN; }

	public TertiaryStructureBondsParser(TokenStream input) {
		super(input);
		_interp = new ParserATNSimulator(this,_ATN,_decisionToDFA,_sharedContextCache);
	}

	public static class BondsContext extends ParserRuleContext {
		public BondsContext(ParserRuleContext parent, int invokingState) {
			super(parent, invokingState);
		}
		@Override public int getRuleIndex() { return RULE_bonds; }
	 
		public BondsContext() { }
		public void copyFrom(BondsContext ctx) {
			super.copyFrom(ctx);
		}
	}
	public static class BondsEndContext extends BondsContext {
		public BondContext bond() {
			return getRuleContext(BondContext.class,0);
		}
		public BondsEndContext(BondsContext ctx) { copyFrom(ctx); }
		@Override
		public void enterRule(ParseTreeListener listener) {
			if ( listener instanceof TertiaryStructureBondsListener ) ((TertiaryStructureBondsListener)listener).enterBondsEnd(this);
		}
		@Override
		public void exitRule(ParseTreeListener listener) {
			if ( listener instanceof TertiaryStructureBondsListener ) ((TertiaryStructureBondsListener)listener).exitBondsEnd(this);
		}
		@Override
		public <T> T accept(ParseTreeVisitor<? extends T> visitor) {
			if ( visitor instanceof TertiaryStructureBondsVisitor ) return ((TertiaryStructureBondsVisitor<? extends T>)visitor).visitBondsEnd(this);
			else return visitor.visitChildren(this);
		}
	}
	public static class BondsContinueContext extends BondsContext {
		public BondContext bond() {
			return getRuleContext(BondContext.class,0);
		}
		public BondsContext bonds() {
			return getRuleContext(BondsContext.class,0);
		}
		public BondsContinueContext(BondsContext ctx) { copyFrom(ctx); }
		@Override
		public void enterRule(ParseTreeListener listener) {
			if ( listener instanceof TertiaryStructureBondsListener ) ((TertiaryStructureBondsListener)listener).enterBondsContinue(this);
		}
		@Override
		public void exitRule(ParseTreeListener listener) {
			if ( listener instanceof TertiaryStructureBondsListener ) ((TertiaryStructureBondsListener)listener).exitBondsContinue(this);
		}
		@Override
		public <T> T accept(ParseTreeVisitor<? extends T> visitor) {
			if ( visitor instanceof TertiaryStructureBondsVisitor ) return ((TertiaryStructureBondsVisitor<? extends T>)visitor).visitBondsContinue(this);
			else return visitor.visitChildren(this);
		}
	}

	public final BondsContext bonds() throws RecognitionException {
		BondsContext _localctx = new BondsContext(_ctx, getState());
		enterRule(_localctx, 0, RULE_bonds);
		try {
			setState(9);
			_errHandler.sync(this);
			switch ( getInterpreter().adaptivePredict(_input,0,_ctx) ) {
			case 1:
				_localctx = new BondsContinueContext(_localctx);
				enterOuterAlt(_localctx, 1);
				{
				setState(4);
				bond();
				setState(5);
				match(T__0);
				setState(6);
				bonds();
				}
				break;
			case 2:
				_localctx = new BondsEndContext(_localctx);
				enterOuterAlt(_localctx, 2);
				{
				setState(8);
				bond();
				}
				break;
			}
		}
		catch (RecognitionException re) {
			_localctx.exception = re;
			_errHandler.reportError(this, re);
			_errHandler.recover(this, re);
		}
		finally {
			exitRule();
		}
		return _localctx;
	}

	public static class BondContext extends ParserRuleContext {
		public List<TerminalNode> INDEX() { return getTokens(TertiaryStructureBondsParser.INDEX); }
		public TerminalNode INDEX(int i) {
			return getToken(TertiaryStructureBondsParser.INDEX, i);
		}
		public BondContext(ParserRuleContext parent, int invokingState) {
			super(parent, invokingState);
		}
		@Override public int getRuleIndex() { return RULE_bond; }
		@Override
		public void enterRule(ParseTreeListener listener) {
			if ( listener instanceof TertiaryStructureBondsListener ) ((TertiaryStructureBondsListener)listener).enterBond(this);
		}
		@Override
		public void exitRule(ParseTreeListener listener) {
			if ( listener instanceof TertiaryStructureBondsListener ) ((TertiaryStructureBondsListener)listener).exitBond(this);
		}
		@Override
		public <T> T accept(ParseTreeVisitor<? extends T> visitor) {
			if ( visitor instanceof TertiaryStructureBondsVisitor ) return ((TertiaryStructureBondsVisitor<? extends T>)visitor).visitBond(this);
			else return visitor.visitChildren(this);
		}
	}

	public final BondContext bond() throws RecognitionException {
		BondContext _localctx = new BondContext(_ctx, getState());
		enterRule(_localctx, 2, RULE_bond);
		try {
			enterOuterAlt(_localctx, 1);
			{
			setState(11);
			match(T__1);
			setState(12);
			match(INDEX);
			setState(13);
			match(T__2);
			setState(14);
			match(INDEX);
			setState(15);
			match(T__3);
			}
		}
		catch (RecognitionException re) {
			_localctx.exception = re;
			_errHandler.reportError(this, re);
			_errHandler.recover(this, re);
		}
		finally {
			exitRule();
		}
		return _localctx;
	}

	public static final String _serializedATN =
		"\u0004\u0001\u0007\u0012\u0002\u0000\u0007\u0000\u0002\u0001\u0007\u0001"+
		"\u0001\u0000\u0001\u0000\u0001\u0000\u0001\u0000\u0001\u0000\u0003\u0000"+
		"\n\b\u0000\u0001\u0001\u0001\u0001\u0001\u0001\u0001\u0001\u0001\u0001"+
		"\u0001\u0001\u0001\u0001\u0000\u0000\u0002\u0000\u0002\u0000\u0000\u0010"+
		"\u0000\t\u0001\u0000\u0000\u0000\u0002\u000b\u0001\u0000\u0000\u0000\u0004"+
		"\u0005\u0003\u0002\u0001\u0000\u0005\u0006\u0005\u0001\u0000\u0000\u0006"+
		"\u0007\u0003\u0000\u0000\u0000\u0007\n\u0001\u0000\u0000\u0000\b\n\u0003"+
		"\u0002\u0001\u0000\t\u0004\u0001\u0000\u0000\u0000\t\b\u0001\u0000\u0000"+
		"\u0000\n\u0001\u0001\u0000\u0000\u0000\u000b\f\u0005\u0002\u0000\u0000"+
		"\f\r\u0005\u0005\u0000\u0000\r\u000e\u0005\u0003\u0000\u0000\u000e\u000f"+
		"\u0005\u0005\u0000\u0000\u000f\u0010\u0005\u0004\u0000\u0000\u0010\u0003"+
		"\u0001\u0000\u0000\u0000\u0001\t";
	public static final ATN _ATN =
		new ATNDeserializer().deserialize(_serializedATN.toCharArray());
	static {
		_decisionToDFA = new DFA[_ATN.getNumberOfDecisions()];
		for (int i = 0; i < _ATN.getNumberOfDecisions(); i++) {
			_decisionToDFA[i] = new DFA(_ATN.getDecisionState(i), i);
		}
	}
}