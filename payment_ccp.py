"""
CCP Payment Integration Module for Streamlit
Integrates Algerian CCP (Compte Chèque Postal) payment functionality
"""

import streamlit as st
import hashlib
import time
from datetime import datetime
import json

# ============================================================================
# PAYMENT CONFIGURATION
# ============================================================================

# CCP Payment Plans
PAYMENT_PLANS = {
    "free": {
        "name": "Free Plan",
        "price": 0,
        "features": [
            "5 map generations per day",
            "Basic satellite imagery",
            "Standard resolution"
        ],
        "map_limit": 5
    },
    "basic": {
        "name": "Basic Plan",
        "price": 500,  # DZD
        "features": [
            "50 map generations per day",
            "All satellite datasets",
            "High resolution imagery",
            "Custom RGB bands",
            "Email support"
        ],
        "map_limit": 50,
        "duration_days": 30
    },
    "pro": {
        "name": "Pro Plan",
        "price": 1500,  # DZD
        "features": [
            "Unlimited map generations",
            "Priority processing",
            "All features included",
            "Bulk downloads",
            "API access",
            "Priority support"
        ],
        "map_limit": -1,  # Unlimited
        "duration_days": 30
    }
}

# CCP Account Information (Replace with your actual CCP details)
CCP_ACCOUNT = {
    "account_number": "0020012345678901",  # 18-digit CCP account
    "account_holder": "Your Name or Business Name",
    "key": "12"  # CCP Key (2 digits)
}

# ============================================================================
# PAYMENT STATE MANAGEMENT
# ============================================================================

def init_payment_state():
    """Initialize payment-related session state"""
    if "user_plan" not in st.session_state:
        st.session_state.user_plan = "free"
    
    if "maps_generated_today" not in st.session_state:
        st.session_state.maps_generated_today = 0
    
    if "last_reset_date" not in st.session_state:
        st.session_state.last_reset_date = datetime.now().date()
    
    if "plan_expiry" not in st.session_state:
        st.session_state.plan_expiry = None
    
    if "pending_payments" not in st.session_state:
        st.session_state.pending_payments = []


def reset_daily_limit():
    """Reset daily map generation counter if new day"""
    today = datetime.now().date()
    if st.session_state.last_reset_date != today:
        st.session_state.maps_generated_today = 0
        st.session_state.last_reset_date = today


def check_usage_limit():
    """Check if user has reached their usage limit"""
    reset_daily_limit()
    
    plan = PAYMENT_PLANS[st.session_state.user_plan]
    limit = plan["map_limit"]
    
    if limit == -1:  # Unlimited
        return True
    
    return st.session_state.maps_generated_today < limit


def increment_usage():
    """Increment map generation counter"""
    st.session_state.maps_generated_today += 1


def get_remaining_uses():
    """Get remaining map generations for today"""
    reset_daily_limit()
    plan = PAYMENT_PLANS[st.session_state.user_plan]
    limit = plan["map_limit"]
    
    if limit == -1:
        return "Unlimited"
    
    remaining = limit - st.session_state.maps_generated_today
    return max(0, remaining)


# ============================================================================
# PAYMENT PROCESSING
# ============================================================================

def generate_payment_reference():
    """Generate unique payment reference number"""
    timestamp = str(int(time.time()))
    random_component = hashlib.md5(timestamp.encode()).hexdigest()[:8]
    return f"SAT-{timestamp[-6:]}-{random_component.upper()}"


def create_payment_record(plan_id, amount, reference):
    """Create a payment record"""
    payment = {
        "reference": reference,
        "plan_id": plan_id,
        "amount": amount,
        "currency": "DZD",
        "status": "pending",
        "created_at": datetime.now().isoformat(),
        "payment_method": "CCP"
    }
    
    st.session_state.pending_payments.append(payment)
    return payment


def verify_payment_manual(reference):
    """Manual payment verification (admin confirms payment received)"""
    # In a real application, this would check against a database
    # For now, we'll use a simple confirmation mechanism
    st.info("⏳ Payment verification pending. Admin will confirm receipt.")
    return False


def activate_plan(plan_id):
    """Activate a paid plan for the user"""
    st.session_state.user_plan = plan_id
    
    if plan_id != "free":
        plan = PAYMENT_PLANS[plan_id]
        duration = plan.get("duration_days", 30)
        from datetime import timedelta
        expiry = datetime.now() + timedelta(days=duration)
        st.session_state.plan_expiry = expiry.isoformat()
    
    st.session_state.maps_generated_today = 0


# ============================================================================
# PAYMENT UI COMPONENTS
# ============================================================================

def display_pricing_cards():
    """Display pricing plans in card format"""
    st.subheader("💳 Choose Your Plan")
    
    cols = st.columns(3)
    
    for idx, (plan_id, plan) in enumerate(PAYMENT_PLANS.items()):
        with cols[idx]:
            # Card styling
            is_current = st.session_state.user_plan == plan_id
            border_color = "#4CAF50" if is_current else "#ddd"
            
            st.markdown(f"""
            <div style="border: 2px solid {border_color}; border-radius: 10px; padding: 20px; height: 100%;">
                <h3 style="color: #1f77b4;">{plan['name']}</h3>
                <h2 style="color: #2ecc71;">{plan['price']:,} DZD</h2>
                <p style="color: #666;">{'per month' if plan_id != 'free' else 'forever'}</p>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("**Features:**")
            for feature in plan["features"]:
                st.markdown(f"✓ {feature}")
            
            if is_current:
                st.success("✅ Current Plan")
            elif plan_id == "free":
                if st.button(f"Downgrade to Free", key=f"btn_{plan_id}", use_container_width=True):
                    activate_plan("free")
                    st.rerun()
            else:
                if st.button(f"Subscribe - {plan['price']:,} DZD", key=f"btn_{plan_id}", 
                           type="primary", use_container_width=True):
                    st.session_state.selected_plan_for_payment = plan_id
                    st.rerun()


def display_ccp_payment_form(plan_id):
    """Display CCP payment instructions and form"""
    plan = PAYMENT_PLANS[plan_id]
    reference = generate_payment_reference()
    
    st.subheader(f"💰 Pay for {plan['name']}")
    st.markdown(f"**Amount to pay:** {plan['price']:,} DZD")
    
    # Display CCP account information
    st.info("📮 **Payment via CCP (Compte Chèque Postal)**")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"""
        **CCP Account Details:**
        - **Account Number:** `{CCP_ACCOUNT['account_number']}`
        - **Key:** `{CCP_ACCOUNT['key']}`
        - **Account Holder:** {CCP_ACCOUNT['account_holder']}
        """)
    
    with col2:
        st.markdown(f"""
        **Your Payment Reference:**
        - **Reference Number:** `{reference}`
        - **Amount:** {plan['price']:,} DZD
        
        ⚠️ **Important:** Include this reference in your payment!
        """)
    
    st.markdown("---")
    
    # Payment instructions
    with st.expander("📋 How to Pay via CCP", expanded=True):
        st.markdown("""
        **Step 1: Go to your nearest Post Office**
        - Bring your CCP card or passbook
        - Or use BaridiMob mobile app
        
        **Step 2: Make the transfer**
        - Transfer to the account number above
        - Include the payment reference in the description
        - Keep your receipt
        
        **Step 3: Submit payment confirmation**
        - Upload your payment receipt below
        - Or enter your CCP transaction reference
        - We'll verify and activate your plan within 24 hours
        """)
    
    st.markdown("---")
    
    # Payment confirmation form
    st.subheader("✅ Confirm Your Payment")
    
    with st.form(f"payment_form_{reference}"):
        payment_method = st.radio(
            "How did you pay?",
            ["CCP Transfer at Post Office", "BaridiMob App", "CCP Online"]
        )
        
        transaction_ref = st.text_input(
            "CCP Transaction Reference",
            placeholder="Enter your CCP transaction reference",
            help="Find this on your receipt or in BaridiMob"
        )
        
        receipt_upload = st.file_uploader(
            "Upload Payment Receipt (Optional)",
            type=["jpg", "jpeg", "png", "pdf"],
            help="Upload a photo or scan of your payment receipt"
        )
        
        sender_ccp = st.text_input(
            "Your CCP Account Number (Optional)",
            placeholder="Your 18-digit CCP account",
            help="This helps us verify the payment faster"
        )
        
        notes = st.text_area(
            "Additional Notes (Optional)",
            placeholder="Any additional information about your payment"
        )
        
        submitted = st.form_submit_button("🚀 Submit Payment Confirmation", type="primary")
        
        if submitted:
            if not transaction_ref:
                st.error("Please enter your CCP transaction reference")
            else:
                # Create payment record
                payment_record = create_payment_record(plan_id, plan['price'], reference)
                payment_record['transaction_ref'] = transaction_ref
                payment_record['payment_method_detail'] = payment_method
                payment_record['sender_ccp'] = sender_ccp
                payment_record['notes'] = notes
                
                # Store receipt if uploaded
                if receipt_upload:
                    payment_record['receipt_uploaded'] = True
                
                st.success("✅ Payment confirmation submitted!")
                st.info("""
                **What's Next?**
                - Our team will verify your payment within 24 hours
                - You'll receive a confirmation once verified
                - Your plan will be activated automatically
                
                **Reference Number:** `{}`
                
                Please keep this reference for your records.
                """.format(reference))
                
                # In a real app, send this to a database and notify admin
                st.balloons()


def display_payment_status():
    """Display current plan status and pending payments"""
    st.subheader("📊 Your Account Status")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        current_plan = PAYMENT_PLANS[st.session_state.user_plan]
        st.metric("Current Plan", current_plan['name'])
    
    with col2:
        remaining = get_remaining_uses()
        st.metric("Maps Remaining Today", remaining)
    
    with col3:
        if st.session_state.plan_expiry:
            expiry = datetime.fromisoformat(st.session_state.plan_expiry)
            days_left = (expiry - datetime.now()).days
            st.metric("Days Until Renewal", days_left)
        else:
            st.metric("Plan Status", "Active")
    
    # Show pending payments
    if st.session_state.pending_payments:
        st.markdown("---")
        st.subheader("⏳ Pending Payments")
        
        for payment in st.session_state.pending_payments:
            if payment['status'] == 'pending':
                with st.expander(f"Payment {payment['reference']}", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**Plan:** {PAYMENT_PLANS[payment['plan_id']]['name']}")
                        st.write(f"**Amount:** {payment['amount']:,} DZD")
                        st.write(f"**Status:** {payment['status'].upper()}")
                    with col2:
                        st.write(f"**Reference:** {payment['reference']}")
                        st.write(f"**Date:** {payment['created_at'][:10]}")
                        st.write(f"**Method:** {payment['payment_method']}")


def show_payment_gate():
    """Show payment gate when usage limit is reached"""
    st.warning("⚠️ You've reached your daily limit!")
    
    st.markdown("""
    ### 🚀 Upgrade to Continue
    
    You've used all your free map generations for today. Upgrade to a paid plan to:
    - Generate more maps immediately
    - Access high-resolution imagery
    - Use advanced features
    - Get priority support
    """)
    
    if st.button("💳 View Plans & Upgrade", type="primary"):
        st.session_state.show_payment_modal = True
        st.rerun()


# ============================================================================
# MAIN PAYMENT INTERFACE
# ============================================================================

def show_payment_interface():
    """Main payment interface"""
    init_payment_state()
    
    st.title("💳 Subscription & Payments")
    
    # Show current status
    display_payment_status()
    
    st.markdown("---")
    
    # Show pricing cards
    display_pricing_cards()
    
    # Show payment form if plan selected
    if hasattr(st.session_state, 'selected_plan_for_payment'):
        st.markdown("---")
        display_ccp_payment_form(st.session_state.selected_plan_for_payment)
        
        if st.button("← Back to Plans"):
            del st.session_state.selected_plan_for_payment
            st.rerun()


# ============================================================================
# HELPER FUNCTIONS FOR MAIN APP INTEGRATION
# ============================================================================

def check_and_enforce_limits():
    """Check usage limits before allowing map generation"""
    init_payment_state()
    
    if not check_usage_limit():
        show_payment_gate()
        return False
    
    return True


def record_map_generation():
    """Record a map generation event"""
    increment_usage()


# ============================================================================
# ADMIN FUNCTIONS (For payment verification)
# ============================================================================

def admin_payment_verification():
    """Admin interface for verifying payments"""
    st.title("🔐 Admin: Payment Verification")
    
    st.warning("⚠️ Admin Access Required")
    
    password = st.text_input("Admin Password", type="password")
    
    if password == "admin123":  # Change this in production!
        st.success("✅ Admin access granted")
        
        if st.session_state.pending_payments:
            for idx, payment in enumerate(st.session_state.pending_payments):
                if payment['status'] == 'pending':
                    with st.expander(f"Payment {payment['reference']}", expanded=True):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.json(payment)
                        
                        with col2:
                            if st.button(f"✅ Approve Payment", key=f"approve_{idx}"):
                                payment['status'] = 'approved'
                                activate_plan(payment['plan_id'])
                                st.success("Payment approved and plan activated!")
                                st.rerun()
                            
                            if st.button(f"❌ Reject Payment", key=f"reject_{idx}"):
                                payment['status'] = 'rejected'
                                st.error("Payment rejected")
                                st.rerun()
        else:
            st.info("No pending payments to verify")
    elif password:
        st.error("❌ Invalid password")
